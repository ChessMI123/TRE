library(dplyr)
library(survival)
library(splines)

adsl <- adsl %>%
    mutate(groups = if_else(group == "controls", 0, 1))

# Prediction grid
times <- seq(min(adsl$CVDyrs), max(adsl$CVDyrs), length.out = 100)

# Define fixed knots based on event times
spline_df <- 2

knot_times <- adsl$CVDyrs[adsl$CVD == 1]

k <- quantile(knot_times,
              probs = seq(0, 1, length.out = spline_df + 1))

k_inner <- k[-c(1, length(k))]
k_bound <- k[c(1, length(k))]

# Cox model with time-varying treatment effect
fit_spline_model <- coxph(
    Surv(time = CVDyrs, event = CVD) ~
        relevel(factor(groups), ref = "0") + tt(groups),
    data = adsl,
    tt = function(x, time, ...) {
        ns(time, knots = k_inner,
           Boundary.knots = k_bound) * x
    }
)

# Extract coefficients and variance-covariance matrix
b <- coef(fit_spline_model)
V <- vcov(fit_spline_model)

# Function to calculate HR and 95% CI
hazard_ratio <- function(t) {

    X <- cbind(1, ns(t, knots = k_inner,
                     Boundary.knots = k_bound))

    logHR <- as.numeric(X %*% b)

    SE_LOGHR <- sqrt(rowSums((X %*% V) * X))

    data.frame(
        time  = t,
        HR    = exp(logHR),
        lower = exp(logHR - qnorm(0.975) * SE_LOGHR),
        upper = exp(logHR + qnorm(0.975) * SE_LOGHR)
    )
}

# HR over time
hr_df <- hazard_ratio(times)
load(file.path(res, "adsl.RData"))

FAS <- adsl %>% filter(FAS == 1)

tab3 <- FAS %>%
  mutate(
    censor_reason = case_when(
      cens == deathdt       ~ "Death",
      cens == treatment_end ~ "End of treatment",
      cens == endofFU       ~ "End of follow-up",
      TRUE                 ~ NA_character_
    ),
    censor_reason = factor(
      censor_reason,
      levels = c("Death", "End of treatment", "End of follow-up")
    )
  )

T3 <- as.data.frame(
  table(tab3$censor_reason, tab3$treatment, useNA = "always")
)

names(T3) <- c("censor_reason", "treatment", "n")

T3$order <- 1:nrow(T3)

write_xlsx(T3, path = file.path(res, "T3.xlsx"))




Pop <- Pop %>%
  mutate(
    censor_reason = case_when(
      cens == deathdt       ~ "Death",
      cens == treatment_end ~ "End of treatment",
      cens == endofFU       ~ "End of follow-up",
      TRUE                 ~ NA_character_
    )
  )

table(Pop$censor_reason, Pop$AC, useNA = "always")





events <- open_dataset(file.path(der, "EVENTS_VTE.parquet")) %>%
  rename(INDATUM_VTE_new = INDATUM_VTE) %>%
  collect()

check <- Pop %>%
  select(LopNr, INDATUM_VTE, baseline) %>%
  inner_join(events, by = "LopNr", relationship = "many-to-many") %>%
  filter(INDATUM_VTE_new > INDATUM_VTE,
         INDATUM_VTE_new <= baseline) %>%
  arrange(LopNr, INDATUM_VTE_new)

View(check)

n_distinct(check$LopNr)



get_disease(Disease = "Arterial_thrombosis", Register = c("SV"), Position = 1, ICD10 = "^(I63|I64|G459|I2[0-5])")
events <- open_dataset(file.path(der, "EVENTS_Arterial_thrombosis.parquet")) %>%
  collect()

check <- Pop %>%
  select(LopNr, INDATUM_VTE, baseline) %>%
  inner_join(events, by = "LopNr", relationship = "many-to-many") %>%
  filter(INDATUM_Arterial_thrombosis >= INDATUM_VTE,
         INDATUM_Arterial_thrombosis <= baseline) %>%
  arrange(LopNr, INDATUM_Arterial_thrombosis)

View(check)

n_distinct(check$LopNr)


time_to_event <- function(data, Name) {

  ind_col    <- paste0("INDATUM_", Name)
  months_col <- paste0(Name, "months")

  # Read event dataset
  events <- open_dataset(file.path(der, paste0("EVENTS_", Name, ".parquet"))) %>%
    select(LopNr, all_of(ind_col)) %>%
    collect()

  # Identify first event after baseline
  first_event <- data %>%
    select(LopNr, baseline) %>%
    inner_join(events, by = "LopNr", relationship = "many-to-many") %>%
    filter(.data[[ind_col]] > baseline) %>%
    group_by(LopNr) %>%
    slice_min(order_by = .data[[ind_col]], n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(LopNr, all_of(ind_col))

  # Calculate time to event in months
  data %>%
    left_join(first_event, by = "LopNr") %>%
    mutate(
      !!Name := case_when(
        !is.na(.data[[ind_col]]) & .data[[ind_col]] <= cens ~ 1,
        TRUE ~ 0
      ),

      !!months_col := case_when(
        .data[[Name]] == 1 ~ as.numeric(.data[[ind_col]] - baseline) / (365.25 / 12),
        TRUE ~ as.numeric(cens - baseline) / (365.25 / 12)
      )
    )
}










png(file.path(projroot,"Output/",paste0("Figure_03_01_",format(Sys.Date(),"%Y-%m-%d"), ".png", sep = "")), width = 7.5, height = 3.5, units="in", res = 300)
p
dev.off()

# Descriptive statistics
desc <- data %>%
  group_by(.data[[group]]) %>%
  summarise(
    N = sum(!is.na(.data[[outcome]])),
    events = sum(.data[[outcome]], na.rm = TRUE),
    median_FU = median(.data[[time]], na.rm = TRUE),
    Q1 = quantile(.data[[time]], 0.25, na.rm = TRUE),
    Q3 = quantile(.data[[time]], 0.75, na.rm = TRUE),
    PY = sum(.data[[time]], na.rm = TRUE),
    IR = 100 * events / PY,
    .groups = "drop"
  )




fmt_py <- function(median, Q1, Q3, PY) {
  sprintf("%.1f (%.1f - %.1f) Sum=%.0f", median, Q1, Q3, PY)
}

result <- data.frame(
  outcome = outcome_label,
  group = levels(data[[group]]),
  events = fmt(desc$events[idx], desc$N[idx]),
  PY = fmt_py(desc$median_FU[idx], desc$Q1[idx],
              desc$Q3[idx], desc$PY[idx]),
  IR = round(desc$IR[idx], 2),
  HR = c("", HR),
  pval = c("", p)
)


cox_table(
  data = adsl,
  group = "treatment",
  outcome = "Major_bleeding",
  adj = c("age", "sex"),
  ref = "LMWH"
)


library(dplyr)
library(survival)
library(broom)

# Format events as n/N (%)
fmt <- function(n, N) {
  sprintf("%d/%d (%.1f%%)", n, N, 100 * n / N)
}

# Cox regression table
cox_table <- function(
    data,
    group,
    outcome,
    adj = NULL,
    ref = "LMWH"
) {

  time <- paste0(outcome, "yrs")

  outcome_label <- attr(data[[outcome]], "label")
  if (is.null(outcome_label)) outcome_label <- outcome

  # Define reference group
  data[[group]] <- relevel(factor(data[[group]]), ref = ref)

  # Descriptive statistics
  desc <- data %>%
    group_by(.data[[group]]) %>%
    summarise(
      N = sum(!is.na(.data[[outcome]])),
      events = sum(.data[[outcome]], na.rm = TRUE),
      PY = sum(.data[[time]], na.rm = TRUE),
      IR = 100 * events / PY,
      .groups = "drop"
    )

  # Helper function for Cox regression
  get_cox <- function(vars) {

    formula <- reformulate(
      vars,
      response = paste0(
        "Surv(", time, ", ", outcome, ")"
      )
    )

    coxph(formula, data = data) %>%
      broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term %in% paste0(group, levels(data[[group]])[-1]))
  }

  # Unadjusted Cox regression
  s <- get_cox(group)

  HR <- fmtci(s$estimate, s$conf.low, s$conf.high)
  p <- pval_func(s$p.value)

  # Combine results
  idx <- match(levels(data[[group]]), desc[[group]])

  result <- data.frame(
    outcome = outcome_label,
    group = levels(data[[group]]),
    events = fmt(desc$events[idx], desc$N[idx]),
    PY = round(desc$PY[idx], 0),
    IR = round(desc$IR[idx], 2),
    HR = c("", HR),
    pval = c("", p)
  )

  # Adjusted Cox regression
  if (!is.null(adj)) {

    s_adj <- get_cox(c(group, adj))

    result$HR_adj <- c(
      "",
      fmtci(s_adj$estimate, s_adj$conf.low, s_adj$conf.high)
    )

    result$pval_adj <- c("", pval_func(s_adj$p.value))
  }

  return(result)
}




cancer_vars <- c(
  "bl_ALL", "bl_Breast_cancer", "bl_CML", "bl_GIST",
  "bl_Kidney_cancer", "bl_Lung_cancer", "bl_pNET",
  "bl_Thyroid_cancer", "bl_Urothelial_cancer"
)

# Patients with multiple cancer types
look <- adsl %>%
  mutate(n_cancer = rowSums(across(all_of(cancer_vars), ~ .x == "Yes"),
                            na.rm = TRUE)) %>%
  filter(n_cancer > 1) %>%
  select(LopNr, n_cancer, all_of(cancer_vars))


# Render categorical variables
render_cat <- function(x, ...) {

  n_total <- sum(!is.na(x))
  values <- if (is.factor(x)) levels(x) else unique(na.omit(x))

  if (all(values %in% c("Yes", "No"))) {
    return(c("Yes" = format_n_pct(sum(x == "Yes", na.rm = TRUE), n_total)))
  }

  if (all(values %in% c("Male", "Female"))) {
    return(c("Male" = format_n_pct(sum(x == "Male", na.rm = TRUE), n_total)))
  }

  c(
    "",
    setNames(
      sapply(values, function(v) {
        format_n_pct(sum(x == v, na.rm = TRUE), n_total)
      }),
      values
    )
  )
}


cox_table(
  data = adsl,
  group = "treatment",
  outcome = "Major_bleeding",
  adj = c("age", "sex"),
  ref = "LMWH"
)




# Render categorical variables
render_cat <- function(x, ...) {

  n_total <- sum(!is.na(x))
  values <- if (is.factor(x)) levels(x) else unique(na.omit(x))

  if (all(values %in% c("Yes", "No"))) {
    values <- "Yes"
  } else if (all(values %in% c("Male", "Female"))) {
    values <- "Male"
  }

  c(
    "" = "",
    setNames(
      sapply(values, function(v) {
        format_n_pct(sum(x == v, na.rm = TRUE), n_total)
      }),
      values
    )
  )
}









VTE2 <- open_dataset(file.path(der, "EVENTS_VTE.parquet")) %>%
  transmute(LopNr, INDATUM_VTE2 = INDATUM_VTE) %>%
  collect()

# Identify additional VTE events between first VTE and baseline
look <- Pop %>%
  inner_join(VTE2, by = "LopNr", relationship = "many-to-many") %>%
  filter(
    INDATUM_VTE2 > INDATUM_VTE,
    INDATUM_VTE2 < baseline
  ) %>%
  arrange(LopNr, INDATUM_VTE2)



time_to_event <- function(data, Name) {

  ind_col   <- paste0("INDATUM_", Name)
  years_col <- paste0(Name, "yrs")

  # Read event dataset
  events <- open_dataset(file.path(der, paste0("EVENTS_", Name, ".parquet"))) %>%
    select(LopNr, all_of(ind_col)) %>%
    collect()

  # Identify first event after baseline
  first_event <- data %>%
    select(LopNr, baseline) %>%
    inner_join(events, by = "LopNr", relationship = "many-to-many") %>%
    filter(.data[[ind_col]] > baseline) %>%
    group_by(LopNr) %>%
    slice_min(order_by = .data[[ind_col]], n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(LopNr, all_of(ind_col))

  # Calculate time to event
  data %>%
    left_join(first_event, by = "LopNr") %>%
    mutate(
      !!Name := case_when(
        !is.na(.data[[ind_col]]) & .data[[ind_col]] <= cens ~ 1,
        TRUE ~ 0
      ),

      !!years_col := case_when(
        .data[[Name]] == 1 ~ as.numeric(.data[[ind_col]] - baseline) / 365.25,
        TRUE ~ as.numeric(cens - baseline) / 365.25
      )
    )
}





time_to_event <- function(data, Name) {

  ind_col   <- paste0("INDATUM_", Name)
  years_col <- paste0(Name, "yrs")

  # Read event dataset
  events <- open_dataset(file.path(der, paste0("EVENTS_", Name, ".parquet"))) %>%
    select(LopNr, all_of(ind_col)) %>%
    collect()

  # Identify first event after baseline
  first_event <- data %>%
    select(LopNr, baseline) %>%
    inner_join(events, by = "LopNr", relationship = "many-to-many") %>%
    filter(.data[[ind_col]] > baseline) %>%
    group_by(LopNr) %>%
    summarise(
      !!ind_col := min(.data[[ind_col]]),
      .groups = "drop"
    )

  # Calculate time to event
  data %>%
    left_join(first_event, by = "LopNr") %>%
    mutate(
      !!Name := case_when(
        !is.na(.data[[ind_col]]) & .data[[ind_col]] <= cens ~ 1,
        TRUE ~ 0
      ),

      !!years_col := case_when(
        .data[[Name]] == 1 ~ as.numeric(.data[[ind_col]] - baseline) / 365.25,
        TRUE ~ as.numeric(cens - baseline) / 365.25
      )
    )
}






Fatal_Bleeding <- Major_bleeding %>%
  filter(REG == "CDR") %>%
  rename(
    Fatal_Bleeding = Major_bleeding,
    INDATUM_Fatal_Bleeding = INDATUM_Major_bleeding
  )

write_parquet(Fatal_Bleeding, file.path(der, "EVENTS_Fatal_Bleeding.parquet"))






render_cat <- function(x, ...) {

  n_total <- sum(!is.na(x))
  values <- if (is.factor(x)) levels(x) else unique(na.omit(x))

  if (all(values %in% c("Yes", "No"))) {
    values <- "Yes"
  } else if (all(values %in% c("Male", "Female"))) {
    values <- "Male"
  }

  setNames(
    sapply(values, function(v) {
      format_n_pct(sum(x == v, na.rm = TRUE), n_total)
    }),
    values
  )
}




c("bl_7y_Hypertension", "bl_7y_HF", "bl_7y_COPD", "bl_7y_Diabetes",
  "bl_7y_Renal", "bl_7y_Liver", "bl_7y_Ischemic_heart_disease",
  "bl_7y_Peripheral", "bl_7y_Ischemicstroke", "bl_7y_Hemorrhagicstroke",
  "bl_7y_CTD", "bl_7y_IBD", "bl_7y_PUD", "bl_7y_AFib")


# Format n (%)
format_n_pct <- function(n, total) {
  sprintf("%d (%.1f%%)", n, 100 * n / total)
}

# Render categorical variables
render_cat <- function(x, ...) {

  n_total <- sum(!is.na(x))

  if (all(na.omit(unique(as.character(x))) %in% c("Yes", "No"))) {
    return(c("Yes" = format_n_pct(sum(x == "Yes", na.rm = TRUE), n_total)))
  }

  if (all(na.omit(unique(as.character(x))) %in% c("Male", "Female"))) {
    return(c("Male" = format_n_pct(sum(x == "Male", na.rm = TRUE), n_total)))
  }

  table1::render.categorical.default(x, ...)
}

c("bl_ALL", "bl_Breast_cancer", "bl_CML", "bl_GIST", "bl_Kidney_cancer",
  "bl_Lung_cancer", "bl_pNET", "bl_Thyroid_cancer", "bl_Urothelial_cancer")
excl_factor = factor(
  excl,
  levels = 0:4,
  labels = c(
    "Included",
    "Dialysis within 3 months or at baseline",
    "Strong CYP3A4 inhibitors within 3 months or at baseline",
    "DOAC or LMWH within 3 months prior to baseline",
    "Same-day initiation of DOAC and LMWH"
  )
)


render_cat <- function(x, ...) {
  if (all(na.omit(unique(as.character(x))) %in% c("Yes", "No"))) {
    n_yes <- sum(x == "Yes", na.rm = TRUE)
    n_total <- sum(!is.na(x))
    
    return(c("Yes" = sprintf("%d (%.1f%%)", n_yes, 100 * n_yes / n_total)))
  }
  
  if (all(na.omit(unique(as.character(x))) %in% c("Male", "Female"))) {
    n_male <- sum(x == "Male", na.rm = TRUE)
    n_total <- sum(!is.na(x))
    
    return(c("Male" = sprintf("%d (%.1f%%)", n_male, 100 * n_male / n_total)))
  }
  
  table1::render.categorical.default(x, ...)
}



### Add treatment vars here!
all_vars <- c("gestational_age_weeks", "birth_weight_g", "sex",
              "Ivh_ae_bin", "nec_ae_bin", "pdamed_ae_bin", "sepsis_ae_bin",
              "type1_rop", "rop_treatment","zone1","zone2","postnatal_age_at_inclusion_weeks",
              "postnatal_age_at_study_stop_weeks", "postmenstrual_age_at_inclusion_weeks",
              "postmenstrual_age_at_study_stop_weeks", "duration_of_intervention_weeks")

t1f <- as.formula(paste("~", paste(all_vars, collapse = "+"), "| randomized_group"))

tab1 <- as.data.frame(
  table1::table1(
    t1f,
    data = adsl_PP,
    render.categorical = render_cat,
    render.continuous = c("median [min-max]"),
    render.missing = NULL,
  )
)






adsl <- adsl %>%
  time_to_event("Major_bleeding")

time_to_event <- function(data, Name) {

  ind_col   <- paste0("INDATUM_", Name)
  event_col <- paste0("EVENTDATE_", Name)
  years_col <- paste0(Name, "yrs")

  # Read event dataset
  events <- open_dataset(file.path(der, paste0("EVENTS_", Name, ".parquet"))) %>%
    select(LopNr, all_of(ind_col)) %>%
    collect()

  # Identify first event after baseline
  first_event <- data %>%
    select(LopNr, baseline) %>%
    inner_join(events, by = "LopNr", relationship = "many-to-many") %>%
    filter(.data[[ind_col]] > baseline) %>%
    group_by(LopNr) %>%
    summarise(
      !!event_col := min(.data[[ind_col]]),
      .groups = "drop"
    )

  # Calculate time to event
  data %>%
    left_join(first_event, by = "LopNr") %>%
    mutate(
      !!Name := case_when(
        !is.na(.data[[event_col]]) & .data[[event_col]] <= cens ~ 1,
        TRUE ~ 0
      ),

      !!years_col := case_when(
        .data[[Name]] == 1 ~ as.numeric(.data[[event_col]] - baseline) / 365.25,
        TRUE ~ as.numeric(cens - baseline) / 365.25
      )
    )
}



# Patients with DOAC or LMWH dispensing within 3 months prior to baseline
look <- Pop %>%
  inner_join(AC, by = "LopNr", relationship = "many-to-many") %>%
  filter(
    EDATUM_AC >= baseline %m-% months(3),
    EDATUM_AC < baseline
  ) %>%
  distinct(LopNr)



Pop %>%
  filter(treatment == "DOAC") %>%
  inner_join(LMWH, by = "LopNr") %>%
  filter(EDATUM_LMWH < baseline) %>%
  distinct(LopNr)



diseases <- c(
  "N_Intra_bleeding",
  "T_Intra_bleeding",
  "Hemopericardium",
  "Gastro_bleeding",
  "Airway_bleeding",
  "Anemia_bleeding",
  "Genital_bleeding",
  "Hematuria",
  "Eye_bleeding"
)

Bleeding <- lapply(diseases, function(disease) {

  open_dataset(file.path(der, paste0("EVENTS_", disease, ".parquet"))) %>%
    select(LopNr, Reg, INDATUM, DIAGNOSE) %>%
    collect()

}) %>%
  bind_rows() %>%
  distinct() %>%
  arrange(LopNr, INDATUM)


write_parquet(Bleeding, file.path(der, "EVENTS_Major_bleeding.parquet"))






bleeding_types <- c(
  "Death_bleeding",
  "N_intracranial_bleeding",
  "T_intracranial_bleeding",
  "Hemopericardium",
  "Gastrointestinalbleeding",
  "Airway_bleeding",
  "Anemia_due_to_bleeding",
  "Genital_bleeding",
  "Hematuria",
  "Eye_bleeding"
)

Major_bleeding <- bind_rows(
  lapply(bleeding_types, function(x) {
    open_dataset(file.path(der, paste0("EVENTS_", x, ".parquet"))) %>%
      collect() %>%
      mutate(bleeding_type = x)
  })
) %>%
  arrange(LopNr, INDATUM)

write_parquet(
  Major_bleeding,
  file.path(der, "EVENTS_Major_bleeding.parquet")
)


read <- open_dataset(file.path(der, "EVENTS_Major_bleeding.parquet")) %>%
  collect() %>%
  mutate(
    bleeding_type = case_when(
      if_any(c(HDIA, DIA1:DIA30), ~ grepl("^I60|^I61|^I62", .x)) ~ "N_intracranial_bleeding",
      if_any(c(HDIA, DIA1:DIA30), ~ grepl("^S064|^S065|^S066", .x)) ~ "T_intracranial_bleeding",
      if_any(c(HDIA, DIA1:DIA30), ~ grepl("^I312", .x)) ~ "Hemopericardium",
      if_any(c(HDIA, DIA1:DIA30), ~ grepl("^I850|^I983|^K250|^K252|^K254|^K256|^K260|^K262|^K264|^K266|^K270|^K272|^K274|^K276|^K280|^K282|^K284|^K286|^K290|^K625|^K920|^K921|^K922", .x)) ~ "Gastrointestinalbleeding",
      if_any(c(HDIA, DIA1:DIA30), ~ grepl("^R04|^R58|^J942|^S271|^H922", .x)) ~ "Airway_bleeding",
      if_any(c(HDIA, DIA1:DIA30), ~ grepl("^D629", .x)) ~ "Anemia_due_to_bleeding",
      if_any(c(HDIA, DIA1:DIA30), ~ grepl("^N920|^N921|^N922|^N923|^N924|^N93|^N950|^N421|^N501A", .x)) ~ "Genital_bleeding",
      if_any(c(HDIA, DIA1:DIA30), ~ grepl("^R31", .x)) ~ "Hematuria",
      if_any(c(HDIA, DIA1:DIA30), ~ grepl("^H113|^H313|^H356|^H431|^H450", .x)) ~ "Eye_bleeding",
      TRUE ~ NA_character_
    )
  )








mutate(treatment_end = pmin(END_TKI, END_AC, switch, na.rm = TRUE))

# Population
Pop <- VTE %>%
  inner_join(TKI, by = "LopNr", relationship = "many-to-many") %>%
  inner_join(AC, by = "LopNr", relationship = "many-to-many") %>%
  filter(
    age >= 18,                              # Age >= 18 years
    START_AC >= INDATUM_VTE,                # AC initiated on or after first VTE
    START_AC >= EDATUM_TKI,                 # No future TKI prescriptions
    START_AC <= EDATUM_TKI + 60,            # Prescriptions within 60 days
    START_AC <= UTDATUM_VTE + 30,           # Within 30 days of discharge
    END_TKI >= START_AC,                    # Ongoing TKI treatment at AC initiation
    UTDATUM_VTE < as.Date("2025-10-01")      # End of inclusion period
  ) %>%
  mutate(baseline = START_AC) %>%
  arrange(LopNr, baseline, desc(EDATUM_TKI)) %>%
  group_by(LopNr) %>%
  slice(1) %>%
  ungroup()
# Population
Pop <- VTE %>%
  inner_join(TKI, by = "LopNr", relationship = "many-to-many") %>%
  inner_join(AC, by = "LopNr", relationship = "many-to-many") %>%
  mutate(baseline = START_AC) %>%
  filter(
    age >= 18,                              # Age >= 18
    START_AC >= START_TKI,                  # AC initiated on or after TKI initiation
    START_AC >= INDATUM_VTE,                # AC initiated on or after VTE diagnosis
    EDATUM_TKI <= baseline,                 # No future TKI prescriptions
    EDATUM_AC == baseline,                  # AC prescription at baseline
    EDATUM_AC <= EDATUM_TKI + 60,           # Prescriptions within 60 days
    baseline <= UTDATUM_VTE + 30,           # Within 30 days of discharge
    END_AC >= baseline,                     # AC still ongoing at baseline
    END_TKI >= baseline,                    # TKI still ongoing at baseline
    UTDATUM_VTE < as.Date("2025-10-01")
  ) %>%
  arrange(LopNr, baseline, desc(EDATUM_TKI)) %>%
  group_by(LopNr) %>%
  slice(1) %>%
  ungroup()



mutate(baseline = pmax(EDATUM_AC, EDATUM_TKI)) %>%
filter(
  age >= 18,
  abs(as.numeric(EDATUM_AC - EDATUM_TKI)) <= 60,
  START_AC >= INDATUM_VTE - 45,
  baseline <= UTDATUM_VTE + 30,
  START_AC >= START_TKI - 10,
  END_AC >= baseline,
  END_TKI >= baseline,
  UTDATUM_VTE < as.Date("2025-10-01")
)



####################################################################################################
# ¤ Load adsl
load(file.path(res,"adsl.RData"))

##########################################################
##########################################################
#### table 0 inclusions/exclusions                    ####
##########################################################
##########################################################

### Individuals with missing diagnosis not excluded
T0 <- as.data.frame(table1::table1(~excl_factor, data = adsl))
## clean up
colnames(T0)<- c("incl_excl", "Overall")
T0<- T0[-2,] #remove unnecessary text. IMPORTANT: DONT RUNT THIS ROW MORE THAN ONCE!!
## save to csv to read in sas when creating report
write_xlsx(x = T0, path = file.path(res,"T0.xlsx"))





#############################################################
#############################################################
#### Table 1.1  Baseline table (FAS)                      ####
#############################################################
#############################################################

T1 <- adsl %>% filter(FAS==1)

### Add treatment vars here!, ""
all_vars<- c("sex","age","age_group_factor","education_factor", "income_quintile_factor","Birth_region","trt_period_factor","bl_MI", "bl_angina", "bl_valvular","bl_card_arrest","bl_AFib","bl_str

t1f<- as.formula(paste("~", paste(all_vars, collapse = "+")))


tab1<- as.data.frame(table1::table1(t1f , data = T1, render.continuous=c("mean (sd) ; median [min-max]")))
tab1$order<- 1:nrow(tab1)
write_xlsx(tab1, path = file.path(res,"T1.xlsx"))






Pop <- Pop %>%
  mutate(
    END_AC_adj = if_else(
      !is.na(deathd) & deathd > END_AC & deathd <= END_AC + 5,
      deathd, END_AC
    ),
    treatment_end = pmin(switch, END_TKI, END_AC_adj, na.rm = TRUE)
  )




# Patients newly included in Pop
Pop_added <- Pop %>%
  anti_join(Pop_old, by = "LopNr")

# Patients no longer included in Pop
Pop_removed <- Pop_old %>%
  anti_join(Pop, by = "LopNr")

# Patients with a changed baseline date
Pop_changed <- Pop %>%
  inner_join(
    Pop_old %>% select(LopNr, baseline_old = baseline),
    by = "LopNr"
  ) %>%
  filter(baseline != baseline_old) %>%
  mutate(diff_days = as.numeric(baseline - baseline_old)) %>%
  arrange(LopNr)

# All affected patients
Pop_affected <- bind_rows(
  Pop_added   %>% mutate(change = "Newly included"),
  Pop_removed %>% mutate(change = "Excluded"),
  Pop_changed %>% mutate(change = "Changed baseline")
) %>%
  arrange(LopNr)

# Summary of changes
data.frame(
  Old_population     = nrow(Pop_old),
  New_population     = nrow(Pop),
  Newly_included     = nrow(Pop_added),
  Excluded           = nrow(Pop_removed),
  Changed_baseline   = nrow(Pop_changed),
  Earlier_baseline   = sum(Pop_changed$diff_days < 0),
  Later_baseline     = sum(Pop_changed$diff_days > 0)
)

# Inspect affected patients
View(Pop_added)
View(Pop_removed)
View(Pop_changed)
View(Pop_affected)




library(lubridate)

treatment_periods_TKI <- function(data) {

  data %>%
    filter(!is.na(EDATUM)) %>%
    arrange(LopNr, EDATUM) %>%
    group_by(LopNr) %>%
    mutate(
      previous_end = lag(EDATUM %m+% months(3)),
      new_episode = if_else(
        is.na(previous_end) | EDATUM > previous_end,
        1L, 0L
      ),
      episode = cumsum(new_episode)
    ) %>%
    group_by(LopNr, episode) %>%
    mutate(
      START_TKI = min(EDATUM),
      END_TKI   = max(EDATUM) %m+% months(3) - 1
    ) %>%
    ungroup()
}




get_disease_all(Disease = "Major_bleeding", Register = c("SV", "CDR"), Position = 2,
                ICD10 = paste0(
                  "^I60|^I61|^I62|",                 # Intracranial bleeding, non-traumatic
                  "^S064|^S065|^S066|",              # Intracranial bleeding, traumatic
                  "^I312|",                          # Hemopericardium
                  "^I850|^I983|^K250|^K252|^K254|^K256|^K260|^K262|^K264|^K266|^K270|^K272|^K274|^K276|^K280|^K282|^K284|^K286|^K290|^K625|^K920|^K921|^K922|", # Gastrointestinal bleeding
                  "^R04|^R58|^J942|^S271|^H922|",    # Airway and other bleeding
                  "^D629|",                          # Anemia due to bleeding
                  "^N920|^N921|^N922|^N923|^N924|^N93|^N950|^N421|^N501A|", # Genital bleeding
                  "^R31|",                           # Hematuria
                  "^H113|^H313|^H356|^H431|^H450"    # Bleeding in the eye
                ))





relocate(HDIA, all_of(paste0("DIA", 1:48)), .after = last_col())
test <- parquet_get_dia_events(
  REG = "CDR",
  Name = "Death_bleeding",
  ICD_string = "^I60",
  position = 1,
  Year_min = 2018,
  Year_max = 2026
)

View(test)
# Identify patients in Pop with prophylactic LMWH prescriptions
Pop_prophylactic <- LMWH %>%
  filter(prophylactic == 1) %>%
  semi_join(Pop, by = "LopNr")

# Unique LopNr
unique(Pop_prophylactic$LopNr)

# Number of patients
n_distinct(Pop_prophylactic$LopNr)

View(Pop_prophylactic)




# Patients with a changed baseline date
Pop_changed <- Pop %>%
  select(LopNr, baseline_old = baseline) %>%
  inner_join(
    Pop_new %>% select(LopNr, baseline_new = baseline),
    by = "LopNr"
  ) %>%
  filter(baseline_old != baseline_new) %>%
  mutate(diff_days = as.numeric(baseline_new - baseline_old)) %>%
  arrange(LopNr)

View(Pop_changed)

# Number of patients with a changed baseline
nrow(Pop_changed)

# Patients newly included in Pop_new
Pop_added <- Pop_new %>%
  anti_join(Pop, by = "LopNr")

# Patients no longer included in Pop_new
Pop_removed <- Pop %>%
  anti_join(Pop_new, by = "LopNr")

# Summary
data.frame(
  Old_population   = nrow(Pop),
  New_population   = nrow(Pop_new),
  Changed_baseline = nrow(Pop_changed),
  Newly_included   = nrow(Pop_added),
  Excluded         = nrow(Pop_removed)
)





# Population
Pop <- VTE %>%
  inner_join(TKI, by = "LopNr", relationship = "many-to-many") %>%
  inner_join(AC, by = "LopNr", relationship = "many-to-many") %>%
  mutate(baseline = pmax(treatment_start, EDATUM_TKI)) %>%
  filter(
    age >= 18,
    treatment_start <= EDATUM_TKI + 60,  # Maximum 60 days after TKI
    treatment_start >= INDATUM_VTE - 45, # VTE no later than 45 days after AC initiation
    baseline <= UTDATUM_VTE + 30,        # Concomitant treatment within 30 days of discharge
    treatment_start >= EDATUM_TKI - 10,  # Maximum 10 days full-dose AC before TKI
    UTDATUM_VTE < as.Date("2025-10-01")
  ) %>%
  arrange(LopNr, baseline, EDATUM_TKI, EDATUM_AC) %>%
  group_by(LopNr) %>%
  slice(1) %>%
  ungroup()

# Save old population
Pop_old <- Pop

# Population
Pop <- VTE %>%
  inner_join(TKI, by = "LopNr", relationship = "many-to-many") %>%
  inner_join(AC, by = "LopNr", relationship = "many-to-many") %>%
  mutate(baseline = pmax(EDATUM_AC, EDATUM_TKI)) %>%
  filter(
    age >= 18,
    EDATUM_AC <= EDATUM_TKI + 60,        # Maximum 60 days after TKI
    EDATUM_AC >= INDATUM_VTE - 45,       # VTE no later than 45 days after AC
    baseline <= UTDATUM_VTE + 30,        # Concomitant treatment within 30 days of discharge
    treatment_start >= EDATUM_TKI - 10,  # Maximum 10 days full-dose AC before TKI
    UTDATUM_VTE < as.Date("2025-10-01")
  ) %>%
  arrange(LopNr, baseline, EDATUM_TKI, EDATUM_AC) %>%
  group_by(LopNr) %>%
  slice(1) %>%
  ungroup()


# Compare old and new population

# Newly included patients
Pop_added <- Pop %>%
  anti_join(Pop_old, by = "LopNr")

# No longer included patients
Pop_removed <- Pop_old %>%
  anti_join(Pop, by = "LopNr")

# Patients with a changed baseline
Pop_changed <- Pop %>%
  inner_join(
    Pop_old %>% select(LopNr, baseline_old = baseline),
    by = "LopNr"
  ) %>%
  filter(baseline != baseline_old) %>%
  mutate(diff_days = as.numeric(baseline - baseline_old))

# Patients with TKI initiation after AC prescription
Pop_TKI_after_AC <- Pop %>%
  filter(EDATUM_TKI > EDATUM_AC) %>%
  mutate(days_AC_before_TKI = as.numeric(EDATUM_TKI - EDATUM_AC))

# All affected patients
Pop_affected <- bind_rows(
  Pop_added %>% mutate(change = "Newly included"),
  Pop_removed %>% mutate(change = "Excluded"),
  Pop_changed %>% mutate(change = "Changed baseline")
) %>%
  arrange(LopNr)

# Summary
data.frame(
  Old_population = nrow(Pop_old),
  New_population = nrow(Pop),
  Newly_included = nrow(Pop_added),
  Excluded = nrow(Pop_removed),
  Changed_baseline = nrow(Pop_changed),
  TKI_after_AC = nrow(Pop_TKI_after_AC)
)

# Inspect affected patients
View(Pop_affected)
View(Pop_TKI_after_AC)
# Check if all LopNr are in Pop_excluded
all(exclLMWH %in% Pop_excluded$LopNr)

# LopNr that are in Pop_excluded
intersect(exclLMWH, Pop_excluded$LopNr)

# LopNr that are NOT in Pop_excluded
setdiff(exclLMWH, Pop_excluded$LopNr)

# Number of LopNr that are NOT in Pop_excluded
length(setdiff(exclLMWH, Pop_excluded$LopNr))



# Add switching date
Pop <- Pop %>%
  left_join(
    Pop %>%
      select(LopNr, baseline, treatment) %>%
      inner_join(AC %>% select(LopNr, EDATUM_AC, treatment_AC = treatment), by = "LopNr") %>%
      filter(EDATUM_AC > baseline, treatment_AC != treatment) %>%
      group_by(LopNr) %>%
      summarise(switch = min(EDATUM_AC), .groups = "drop"),
    by = "LopNr"
  )



prophylactic = case_when(
  ATC == "B01AB04" & ANTAL <= 5 ~ 1,
  TRUE ~ NA_real_
)


# Compare baseline dates between Pop and Pop_new
check <- Pop %>%
  select(LopNr, baseline_old = baseline) %>%
  inner_join(Pop_new %>% select(LopNr, baseline_new = baseline), by = "LopNr") %>%
  filter(baseline_old != baseline_new) %>%
  mutate(days_diff = as.numeric(baseline_new - baseline_old)) %>%
  arrange(LopNr)

# Newly included patients
new_patients <- Pop_new %>%
  anti_join(Pop, by = "LopNr")

# Patients no longer included
removed_patients <- Pop %>%
  anti_join(Pop_new, by = "LopNr")

# Summary
nrow(Pop)               # Original population
nrow(Pop_new)           # Updated population
nrow(check)             # Patients with changed baseline
nrow(new_patients)      # Newly included patients
nrow(removed_patients)  # Patients no longer included

# View results
View(check)
View(new_patients)
View(removed_patients)







PopDOAC <- Pop3 %>%
  filter(excl == 0 & treatment == "DOAC")

check <- PopDOAC %>%
  inner_join(LMWH, by = "LopNr") %>%
  arrange(LopNr, EDATUM)


check <- PopDOAC %>%
  inner_join(LMWH, by = "LopNr") %>%
  filter(EDATUM > baseline) %>%
  mutate(days_after_baseline = as.numeric(EDATUM - baseline)) %>%
  arrange(LopNr, EDATUM)


excl3 <- excl3 %>%
  mutate(baseline = if_else(treatment_start <= EDATUM_TKI,
                            pmax(INDATUM_VTE, EDATUM_TKI),
                            baseline))



check <- adsl %>%
  mutate(n_cancer = rowSums(across(c(bl_ALL, bl_Breast_cancer, bl_CML, bl_GIST,
                                    bl_Kidney_cancer, bl_Lung_cancer, bl_pNET,
                                    bl_Thyroid_cancer, bl_Urothelial_cancer),
                                  ~ .x == "Yes"), na.rm = TRUE)) %>%
  filter(n_cancer > 1)


# Demographics
attr(adsl$age, "label") <- "Age at baseline"
attr(adsl$age_group_factor, "label") <- "Age group at baseline"
attr(adsl$Birth_region, "label") <- "Birth region"
attr(adsl$sex, "label") <- "Sex"
attr(adsl$education_factor, "label") <- "Education level"
attr(adsl$income_quintile_factor, "label") <- "Income quintile (Based on KPI adjusted incomes)"
attr(adsl$trt_period_factor, "label") <- "Treatment period"

# Comorbidities (within 7 years prior to or at baseline)
attr(adsl$bl_7y_Hypertension, "label")           <- "Hypertension"
attr(adsl$bl_7y_HF, "label")                     <- "Heart failure"
attr(adsl$bl_7y_COPD, "label")                   <- "Chronic obstructive pulmonary disease"
attr(adsl$bl_7y_Diabetes, "label")               <- "Diabetes mellitus"
attr(adsl$bl_7y_Renal, "label")                  <- "Renal disease"
attr(adsl$bl_7y_Liver, "label")                  <- "Liver disease"
attr(adsl$bl_7y_Ischemic_heart_disease, "label") <- "Ischemic heart disease"
attr(adsl$bl_7y_Peripheral, "label")             <- "Peripheral arterial disease"
attr(adsl$bl_7y_Ischemicstroke, "label")         <- "Ischemic stroke"
attr(adsl$bl_7y_Hemorrhagicstroke, "label")      <- "Hemorrhagic stroke"
attr(adsl$bl_7y_CTD, "label")                    <- "Systemic connective tissue disorders"
attr(adsl$bl_7y_IBD, "label")                    <- "Inflammatory bowel disease"
attr(adsl$bl_7y_PUD, "label")                    <- "Peptic ulcer disease"
attr(adsl$bl_7y_AFib, "label")                   <- "Atrial fibrillation"

# Cancer types (prior to or at baseline)
attr(adsl$bl_ALL, "label")                       <- "Acute lymphatic leukemia (ALL)"
attr(adsl$bl_Breast_cancer, "label")             <- "Breast cancer"
attr(adsl$bl_CML, "label")                       <- "Chronic myeloid leukemia (CML)"
attr(adsl$bl_GIST, "label")                      <- "Gastrointestinal stromal tumor (GIST)"
attr(adsl$bl_Kidney_cancer, "label")             <- "Kidney cancer"
attr(adsl$bl_Lung_cancer, "label")               <- "Lung cancer"
attr(adsl$bl_pNET, "label")                      <- "Pancreatic neuroendocrine tumor (pNET)"
attr(adsl$bl_Thyroid_cancer, "label")            <- "Thyroid cancer"
attr(adsl$bl_Urothelial_cancer, "label")         <- "Urothelial cancer"

# Prior medications (within 6 months prior to or at baseline)
attr(adsl$bl_6m_Antiplatelet, "label")           <- "Antiplatelet treatment"
attr(adsl$bl_6m_PPI, "label")                    <- "Proton pump inhibitors"
attr(adsl$bl_6m_Statins, "label")                <- "Statins"

# Temporary provoking factors (within 3 months prior to or at baseline)
attr(adsl$bl_3m_Nervous_system_surgery, "label")       <- "Surgery on the nervous system including cranial nerves"
attr(adsl$bl_3m_ENT_surgery, "label")                  <- "Ear, nose and throat surgery"
attr(adsl$bl_3m_Heart_vessel_surgery, "label")         <- "Surgery on the heart and large vessels"
attr(adsl$bl_3m_Thoracic_lung_surgery, "label")         <- "Thoracic and lung procedures"
attr(adsl$bl_3m_Gastrointestinal_surgery, "label")      <- "Gastrointestinal surgery"
attr(adsl$bl_3m_Urinary_genital_surgery, "label")       <- "Urinary tract and genital surgery"
attr(adsl$bl_3m_Musculoskeletal_surgery, "label")       <- "Musculoskeletal surgery"
attr(adsl$bl_3m_Obstetrical_surgery, "label")            <- "Obstetrical surgery"
attr(adsl$bl_3m_Central_venous_catheter, "label")       <- "Central venous catheter and implantable venous port"
attr(adsl$bl_3m_Lower_extremity_fracture, "label")      <- "Lower extremity fracture"
attr(adsl$bl_3m_Trauma, "label")                       <- "Trauma"
attr(adsl$bl_3m_Hormone_replacement_therapy, "label")   <- "Hormone replacement therapy"


save(adsl, file = file.path(res,"adsl.RData"))
write_xlsx(adsl, path = file.path(res,"adsl.xlsx"))









# Fix formats and wrap up the data
adsl <- adsl10 %>%
  mutate(sex=factor(Kon, level=1:2, labels=c("Male", "Female")),
         Birth_region = factor(Birth_region, levels = c("Sweden","Nordic (excl Sweden)","Europe (excl Nordic)","Other")),
         excl_factor= factor(excl,
                             levels=0:6,
                             labels=c("Included",
                                      "Excluded due to baseline date not identified",
                                      "Excluded due to heart failure before 2009 or on last day of follow-up",
                                      "Excluded due to death on the same day as heart failure",
                                      "Excluded due to reused personal identification number",
                                      "Excluded due to CIED implantation before baseline",
                                      "Excluded due to first ever recorded CIED after baseline is Utbyte")
         ),
         education_gr= case_when(education_level==1~1,
                                 education_level<5 ~2,
                                 education_level>=5~3,
                                 TRUE~NA),
         education_factor= factor(education_gr, levels=1:3,labels=c("0-9 Years","9-12 Years","More than 12 Years")),

         age_group = case_when(
                     age < 65 ~1,
                     age < 75 ~2,
                     age < 85 ~3,
                     age >= 85 ~4
         ),

         age_group_factor = factor(age_group, levels = 1:4, labels = c("<65","65-74","75-84","¤{unicode '2265'x}85")),

         HF_type = case_when(
           LVEF_PERCENT < 40 ~1,
           LVEF_PERCENT >=40 & LVEF_PERCENT < 50 ~2,
           LVEF_PERCENT >=50  ~3,
           TRUE ~ NA),

         HF_factor = factor(HF_type, levels = 1:3, labels = c("HFrEF","HFmrEF","HFpEF")),

         trt_period = case_when(lisayr <=2012 ~1,
                                lisayr >=2013 & lisayr<=2016~2,
                                lisayr >=2017 & lisayr<=2020~3,
                                lisayr >=2021 & lisayr<=2023~4),

         trt_period_factor= factor(trt_period, levels=1:4, labels=c("2009-2012", "2013-2016", "2017-2020", "2021-2023")),


         income_quintile_factor= factor(income_quintile, levels=1:5,labels=c("1st quintile", "2nd quintile", "3rd quintile", "4th quintile", "5th quintile"))


  )







check <- Pop3 %>%
  filter(treatment == "DOAC", !LopNr %in% excl3$LopNr) %>%
  inner_join(DOAC %>% select(LopNr, EDATUM), by = "LopNr") %>%
  filter(EDATUM >= baseline %m-% months(3), EDATUM < baseline) %>%
  arrange(LopNr, EDATUM)


check <- Pop3 %>%
  filter(treatment == "DOAC", !LopNr %in% excl3$LopNr) %>%
  inner_join(DOAC %>% select(LopNr, EDATUM), by = "LopNr") %>%
  filter(EDATUM >= baseline %m-% months(3), EDATUM < baseline) %>%
  arrange(LopNr, EDATUM)


Pop3 <- Pop2 %>%
  mutate(excl = case_when(excl == 0 & LopNr %in% c(
    1564728, 2499503, 3124008, 3330799, 3338217, 3392967, 3529066, 3544843, # excl DOAC
    110126, 187371, 233004, 504389, 669175, 742191, 774767, 844354,          # excl LMWH
    1721193, 1805284, 1979265, 2383936, 2504734, 2650247, 2680927, 2789045,
    2932662, 2982578, 3009904, 3013373, 3025960, 3455492, 3491378, 3507612,
    3590077, 3616830, 3726200, 3736811, 3830285
  ) ~ 3, TRUE ~ excl)) # We have classified these manually based on their prescription records.





### Temporary provoking factors
adsl6 <- adsl5 %>%
  baseline_3m("Nervous_system_surgery") %>%
  baseline_3m("ENT_surgery") %>%
  baseline_3m("Heart_vessel_surgery") %>%
  baseline_3m("Thoracic_lung_surgery") %>%
  baseline_3m("Gastrointestinal_surgery") %>%
  baseline_3m("Urinary_genital_surgery") %>%
  baseline_3m("Musculoskeletal_surgery") %>%
  baseline_3m("Obstetrical_surgery") %>%
  baseline_3m("Central_venous_catheter") %>%
  baseline_3m("Lower_extremity_fracture") %>%
  baseline_3m("Trauma") %>%
  baseline_3m("Hormone_replacement_therapy")



# Temporary provoking factors - Surgical interventions
get_intervention_all(NAME = "Nervous_system_surgery", KVÅ_CODE = "^AA[A-F][0-9]{2}|^AAH[1-9][0-9]|^AAJ[0-9]{2}|^AAK10|^AAK40|^AAK99|^AA[LMNP][0-9]{2}|^AAU00|^AAW(0[1-9]|[1-9][0-9])")
get_intervention_all(NAME = "ENT_surgery", KVÅ_CODE = "^D[A-HJ][A-W]|^DKW|^DL[D-W]|^D[MNPQW][A-W]")
get_intervention_all(NAME = "Heart_vessel_surgery", KVÅ_CODE = "^F")
get_intervention_all(NAME = "Thoracic_lung_surgery", KVÅ_CODE = "^G[ABDEW][A-W]|^GCA")
get_intervention_all(NAME = "Gastrointestinal_surgery", KVÅ_CODE = "^J")
get_intervention_all(NAME = "Urinary_genital_surgery", KVÅ_CODE = "^K|^L")
get_intervention_all(NAME = "Musculoskeletal_surgery", KVÅ_CODE = "^N")
get_intervention_all(NAME = "Obstetrical_surgery", KVÅ_CODE = "^M")
get_intervention_all(NAME = "Central_venous_catheter", KVÅ_CODE = "^TPH05|^TPH15|^TPH20|^TPX10|^TPW99|^SP(10[5-9]|1[1-9][0-9]|[2-9][0-9]{2})")

# Temporary provoking factors - Diagnoses
get_disease_all(Disease = "Lower_extremity_fracture", ICD10 = "^S72|^S82|^S92|^T12")
get_disease_all(Disease = "Trauma", ICD10 = "^V(0[1-9]|[1-9][0-9])|^W[0-9]{2}|^X[0-5][0-9]")

# Temporary provoking factors - Medications
get_medication_all(NAME = "Hormone_replacement_therapy", ICD_CODE = "G03C|G03F")




### Prior medications
adsl5 <- adsl4 %>%
  baseline_6m("Antiplatelet") %>%
  baseline_6m("PPI") %>%
  baseline_6m("Statins")


get_medication_all(NAME = "Antiplatelet", ICD_CODE = "B01AC")
get_medication_all(NAME = "PPI",          ICD_CODE = "A02BC")
get_medication_all(NAME = "Statins",      ICD_CODE = "C10AA")



baseline <- function(data, Name) {

  message(Name)
  file_path <- file.path(der, paste0("EVENTS_", Name, ".parquet"))

  if (!file.exists(file_path)) {
    file_path <- file.path(meds, paste0("EVENTS_", Name, "_prescription.parquet"))
  }

  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(data)
  }

  tmp <- read_parquet(file_path)

  date_var <- if ("UTDATUM" %in% names(tmp)) "UTDATUM" else if ("EDATUM" %in% names(tmp)) "EDATUM" else paste0("EDATUM_", Name)

  tmp <- tmp %>%
    transmute(LopNr, UTDATUM = .data[[date_var]])

  col_name <- paste0("bl_", Name)
  if (col_name %in% colnames(data)) {data <- data %>% select(-all_of(col_name))}

  df_adsl <- data %>%
    left_join(tmp, by = "LopNr")

  df_temporary <- df_adsl %>%
    filter(UTDATUM <= baseline)

  df_for_join <- data.frame(LopNr = unique(df_temporary$LopNr)) %>%
    mutate(!!sym(col_name) := 1)

  data <- data %>% left_join(df_for_join, by = "LopNr")

  data[col_name][is.na(data[col_name])] <- 0
  data[[col_name]] <- ynf(data[[col_name]])

  return(data)
}



baseline_6m <- function(data, Name) {

  message(Name)
  file_path <- file.path(der, paste0("EVENTS_", Name, ".parquet"))

  if (!file.exists(file_path)) {
    file_path <- file.path(meds, paste0("EVENTS_", Name, "_prescription.parquet"))
  }

  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(data)
  }

  tmp <- read_parquet(file_path)

  date_var <- if ("UTDATUM" %in% names(tmp)) "UTDATUM" else if ("EDATUM" %in% names(tmp)) "EDATUM" else paste0("EDATUM_", Name)

  tmp <- tmp %>%
    transmute(LopNr, UTDATUM = .data[[date_var]])

  col_name <- paste0("bl_6m_", Name)
  if (col_name %in% colnames(data)) {data <- data %>% select(-all_of(col_name))}

  df_adsl <- data %>%
    left_join(tmp, by = "LopNr")

  df_temporary <- df_adsl %>%
    filter(UTDATUM >= baseline %m-% months(6), UTDATUM <= baseline)

  df_for_join <- data.frame(LopNr = unique(df_temporary$LopNr)) %>%
    mutate(!!sym(col_name) := 1)

  data <- data %>% left_join(df_for_join, by = "LopNr")

  data[col_name][is.na(data[col_name])] <- 0
  data[[col_name]] <- ynf(data[[col_name]])

  return(data)
}





### Cancer types
adsl4 <- adsl3 %>%
  baseline("ALL") %>%
  baseline("Breast_cancer") %>%
  baseline("CML") %>%
  baseline("GIST") %>%
  baseline("Kidney_cancer") %>%
  baseline("Lung_cancer") %>%
  baseline("pNET") %>%
  baseline("Thyroid_cancer") %>%
  baseline("Urothelial_cancer")





### Comorbidities
### Comorbidities
adsl3 <- adsl2 %>%
  baseline_7y("Hypertension") %>%
  baseline_7y("HF") %>%
  baseline_7y("COPD") %>%
  baseline_7y("Diabetes") %>%
  baseline_7y("Renal") %>%
  baseline_7y("Liver") %>%
  baseline_7y("Ischemic_heart_disease") %>%
  baseline_7y("Peripheral") %>%
  baseline_7y("Ischemicstroke") %>%
  baseline_7y("Hemorrhagicstroke") %>%
  baseline_7y("CTD") %>%
  baseline_7y("IBD") %>%
  baseline_7y("PUD") %>%
  baseline_7y("AFib")


### Prior Medications
adsl4 <- adsl3 %>%
  baseline_medication("BETA_BLOCKER") %>%
  baseline_medication("RAS_inhibitors") %>%
  baseline_medication("MRA") %>%
  baseline_medication("ARNI") %>%
  baseline_medication("SGLT2")



baseline_7y <- function(data, Name) {

  message(Name)
  file_path <- file.path(der, paste0("EVENTS_", Name, ".parquet"))

  if (!file.exists(file_path)) {
    file_path <- file.path(meds, paste0("EVENTS_", Name, "_prescription.parquet"))
  }

  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(data)
  }

  tmp <- read_parquet(file_path)

  date_var <- if ("UTDATUM" %in% names(tmp)) "UTDATUM" else if ("EDATUM" %in% names(tmp)) "EDATUM" else paste0("EDATUM_", Name)

  tmp <- tmp %>%
    transmute(LopNr, UTDATUM = .data[[date_var]])

  col_name <- paste0("bl_7y_", Name)
  if (col_name %in% colnames(data)) {data <- data %>% select(-all_of(col_name))}

  df_adsl <- data %>%
    left_join(tmp, by = "LopNr")

  df_temporary <- df_adsl %>%
    filter(UTDATUM >= baseline %m-% years(7), UTDATUM <= baseline)

  df_for_join <- data.frame(LopNr = unique(df_temporary$LopNr)) %>%
    mutate(!!sym(col_name) := 1)

  data <- data %>% left_join(df_for_join, by = "LopNr")

  data[col_name][is.na(data[col_name])] <- 0
  data[[col_name]] <- ynf(data[[col_name]])

  return(data)
}





















### Add death
adsl5 <- adsl4 %>%
                  left_join(open_dataset(file.path(der,"CV_death_ULORSAK.parquet"))%>% collect() %>% select(LopNr,CVdeath), by="LopNr") %>%
                  left_join(open_dataset(file.path(der,"SCD_death_ULORSAK.parquet"))%>% collect() %>% select(LopNr,SCDdeath), by="LopNr")


adsl6 <- adsl5 %>% mutate(deathdt = if_else(deathdt <= baseline, baseline, deathdt),
                          endofFU=ymd("2023-12-31"),
                          cens=pmin(endofFU,deathdt, na.rm=TRUE),
                          censyrs = as.numeric(cens - baseline)/365.25,
                          death = case_when(!is.na(deathdt) & deathdt<=endofFU~1,
                                            !is.na(deathdt) & deathdt>endofFU ~0,
                                             is.na(deathdt)                   ~0),
                          deathyrs=as.numeric(case_when(death==1~(deathdt-baseline)/365.25,
                                                        death==0~(endofFU-baseline)/365.25)),
                          deathyrs = if_else(deathyrs == 0,1/365.25,deathyrs),
                          cv_deathyrs_comprisk = deathyrs,
                          cv_death_comprisk = case_when(death == 0 ~ 0,
                                               death == 1 & CVdeath == 1~1,
                                               death == 1 ~ 2),
                          scd_deathyrs_comprisk = deathyrs,
                          scd_death_comprisk = case_when(death == 0 ~ 0,
                                               death == 1 & SCDdeath == 1~1,
                                               death == 1 ~ 2))



Pop3 <- Pop2 %>%
  mutate(excl = case_when(excl == 0 & LopNr %in% c(1564728, 2499503, 3124008, 3330799, 3338217, 3392967, 3529066, 3544843) ~ 3, TRUE ~ excl))

# DOAC
DOACkoll <- DOAC %>%
  filter(LopNr %in% excl3$LopNr) %>%
  left_join(Pop %>% select(LopNr, baseline), by = "LopNr")

# LMWH
LMWHkoll <- LMWH %>%
  filter(LopNr %in% excl3$LopNr) %>%
  left_join(Pop %>% select(LopNr, baseline), by = "LopNr")
baseline_3m <- function(data, Name) {

  message(Name)
  file_path <- file.path(der, paste0("EVENTS_", Name, ".parquet"))

  if (!file.exists(file_path)) {
    file_path <- file.path(meds, paste0("EVENTS_", Name, "_prescription.parquet"))
  }

  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(data)
  }

  tmp <- read_parquet(file_path)

  date_var <- if ("UTDATUM" %in% names(tmp)) "UTDATUM" else if ("EDATUM" %in% names(tmp)) "EDATUM" else paste0("EDATUM_", Name)

  tmp <- tmp %>%
    transmute(LopNr, UTDATUM = .data[[date_var]])

  col_name <- paste0("bl_3m_", Name)
  if (col_name %in% colnames(data)) {data <- data %>% select(-all_of(col_name))}

  df_adsl <- data %>%
    left_join(tmp, by = "LopNr")

  df_temporary <- df_adsl %>%
    filter(UTDATUM >= baseline %m-% months(3), UTDATUM <= baseline)

  df_for_join <- data.frame(LopNr = unique(df_temporary$LopNr)) %>%
    mutate(!!sym(col_name) := 1)

  data <- data %>% left_join(df_for_join, by = "LopNr")

  data[col_name][is.na(data[col_name])] <- 0
  data[[col_name]] <- ynf(data[[col_name]])

  return(data)
}




# Combine all CYP3A4 inhibitor prescriptions
CYP3A4 <- bind_rows(
  open_dataset(file.path(meds, "EVENTS_Itraconazole_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Ketoconazole_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Posaconazole_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Voriconazole_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Ritonavir_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Clarithromycin_prescription.parquet")) %>% collect()
) %>%
  arrange(LopNr, EDATUM) %>%
  transmute(LopNr, EDATUM_CYP3A4 = EDATUM)

write_parquet(CYP3A4, file.path(meds, "EVENTS_CYP3A4_prescription.parquet"))




# CYP3A4 INHIBITORS
get_medication_all(NAME = "Itraconazole",  ICD_CODE = "J02AC02")
get_medication_all(NAME = "Ketoconazole",  ICD_CODE = "J02AB02")
get_medication_all(NAME = "Posaconazole",  ICD_CODE = "J02AC04")
get_medication_all(NAME = "Voriconazole",  ICD_CODE = "J02AC03")
get_medication_all(NAME = "Ritonavir",     ICD_CODE = "J05AE03")
get_medication_all(NAME = "Clarithromycin", ICD_CODE = "J01FA09")



Pop1 <- Pop %>%
  baseline_6m("Dialysis") %>%
  mutate(excl1 = if_else(bl_Dialysis == "Yes", 0, 1))
# Exclusion 1: Dialysis at baseline or within 3 months prior
Pop1 <- Pop %>%
  baseline_disease("Dialysis", uppbound = 0, lowbound = -92) %>%
  mutate(excl1 = if_else(bl_Dialysis == "Yes", 0, 1))




baseline_6m <- function(data, Name) {

  message(Name)
  file_path <- file.path(der, paste0("EVENTS_", Name, ".parquet"))

  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(data)
  }

  tmp <- read_parquet(file_path) %>% select(LopNr, UTDATUM)
  col_name <- paste0("bl_", Name)
  if (col_name %in% colnames(data)) {data <- data %>% select(-all_of(col_name))}

  df_adsl <- data %>%
    left_join(tmp, by = "LopNr")

  df_temporary <- df_adsl %>%
    filter(UTDATUM >= baseline %m-% months(6), UTDATUM <= baseline)

  df_for_join <- data.frame(LopNr = unique(df_temporary$LopNr)) %>%
    mutate(!!sym(col_name) := 1)

  data <- data %>% left_join(df_for_join, by = "LopNr")

  data[col_name][is.na(data[col_name])] <- 0
  data[[col_name]] <- ynf(data[[col_name]])

  return(data)
}

baseline_3m <- function(data, Name) {

  message(Name)
  file_path <- file.path(der, paste0("EVENTS_", Name, ".parquet"))

  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(data)
  }

  tmp <- read_parquet(file_path) %>% select(LopNr, UTDATUM)
  col_name <- paste0("bl_", Name)
  if (col_name %in% colnames(data)) {data <- data %>% select(-all_of(col_name))}

  df_adsl <- data %>%
    left_join(tmp, by = "LopNr")

  df_temporary <- df_adsl %>%
    filter(UTDATUM >= baseline %m-% months(3), UTDATUM <= baseline)

  df_for_join <- data.frame(LopNr = unique(df_temporary$LopNr)) %>%
    mutate(!!sym(col_name) := 1)

  data <- data %>% left_join(df_for_join, by = "LopNr")

  data[col_name][is.na(data[col_name])] <- 0
  data[[col_name]] <- ynf(data[[col_name]])

  return(data)
}

baseline_disease <- function(data, other_diseases, uppbound = 0, lowbound = -Inf) {

  message(other_diseases)
  file_path <- file.path(der, paste0("EVENTS_", gsub("2", "1996", other_diseases), ".parquet"))

  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(data)
  }

  tmp <- read_parquet(file_path) %>% select(LopNr, INDATUM)
  col_name <- paste0("bl_", other_diseases)
  if (col_name %in% colnames(data)) {data <- data %>% select(-all_of(col_name))}

  df_adsl <- data %>%
    left_join(tmp, by = "LopNr")

  df_adsl$days <- as.numeric(as.Date(df_adsl$INDATUM) - as.Date(df_adsl$baseline))

  df_temporary <- df_adsl %>% filter(days <= uppbound, days >= lowbound)
  df_for_join <- data.frame(LopNr = unique(df_temporary$LopNr)) %>%
    mutate(!!sym(col_name) := 1)

  data <- data %>% left_join(df_for_join, by = "LopNr")

  data[col_name][is.na(data[col_name])] <- 0
  data[[col_name]] <- ynf(data[[col_name]])

  return(data)
}





# Exclusion 1: Dialysis at baseline or within 3 months prior
check <- open_dataset(file.path(der, "EVENTS_Dialysis.parquet")) %>%
  collect() %>%
  inner_join(baseline_dat %>% select(LopNr, baseline), by = "LopNr") %>%
  filter(UTDATUM >= baseline %m-% months(3), UTDATUM <= baseline) %>%
  distinct(LopNr)

baseline_dat <- baseline_dat %>%
  mutate(excl = if_else(LopNr %in% check$LopNr, 0, 1))



# Calculate duration of each prescription
DOAC <- DOAC %>%
  mutate(
    tablets = ANTAL * antnum,
    days_supply = case_when(
      ATC == "B01AF01" & styrknum == 15 ~ if_else(tablets <= 42, tablets / 2, 21 + (tablets - 42)), # Rivaroxaban 15 mg: 2 tablets/day for 21 days, thereafter 1/day
      ATC == "B01AF01" & styrknum == 20 ~ tablets, # Rivaroxaban 20 mg: 1 tablet/day
      ATC == "B01AF02" & styrknum == 5 ~ if_else(tablets <= 28, tablets / 4, 7 + (tablets - 28) / 2), # Apixaban 5 mg: 4 tablets/day for 7 days, thereafter 2/day
      ATC == "B01AF02" & styrknum == 2.5 ~ tablets / 2, # Apixaban 2.5 mg: 2 tablets/day
      ATC %in% c("B01AF03", "B01AE07") ~ tablets / 2, # Edoxaban and Dabigatran: 2 tablets/day
      TRUE ~ NA_real_
    ),
    days_supply = ceiling(days_supply),
    supply_end = EDATUM + days_supply - 1
  )






# Combine all DOAC prescriptions
DOAC <- bind_rows(
  open_dataset(file.path(meds, "EVENTS_Rivaroxaban_prescription.parquet")) %>%
    filter(!(coalesce(styrknum %in% c(2.5, 1), FALSE))) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Apixaban_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Dabigatran_prescription.parquet")) %>%
    filter(!(coalesce(styrknum == 75, FALSE))) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Edoxaban_prescription.parquet")) %>% collect()
) %>%
  arrange(LopNr, EDATUM)



adsl %>%
  filter(EDATUM_TKI < INDATUM_VTE) %>%
  mutate(days_before = as.numeric(INDATUM_VTE - EDATUM_TKI)) %>%
  slice_max(days_before, n = 1)



# Combine all AC prescriptions
AC <- bind_rows(
  DOAC %>% transmute(LopNr, AC_date = DOAC_date, treatment = "DOAC"),
  LMWH %>% transmute(LopNr, AC_date = LMWH_date, treatment = "LMWH")
) %>%
  arrange(LopNr, AC_date)

# Save as Parquet
write_parquet(AC, file.path(meds, "EVENTS_AC_prescription.parquet"))



baseline_dat <- VTE %>%
  inner_join(TKI, by = "LopNr", relationship = "one-to-many") %>%
  arrange(LopNr, abs(as.numeric(EDATUM_TKI - INDATUM_VTE))) %>%
  group_by(LopNr) %>%
  slice(1) %>%
  ungroup()


# Combine all TKI prescriptions
TKI <- bind_rows(
  open_dataset(file.path(meds, "EVENTS_Crizotinib_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Ceritinib_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Erdafitinib_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Imatinib_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Osimertinib_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Selpercatinib_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Sotorasib_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Sunitinib_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Tepotinib_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Tucatinib_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Vandetanib_prescription.parquet")) %>% collect()
) %>%
  arrange(LopNr, EDATUM)

filter(INDATUM >= as.Date("2018-01-01") & INDATUM < as.Date("2025-10-01"))
library(arrow)

SV <- open_dataset("/safe/data/Original data 2026/PARQUET/ut_r_par_sv_34759_2025.parquet")

names(SV)

# Calculate duration of each prescription
LMWH <- LMWH %>%
  mutate(
    doses = ANTAL * antnum,
    days_supply = case_when(
      ATC == "B01AB10" ~ doses, # Tinzaparin: 1 dose/day
      ATC == "B01AB04" ~ doses, # Dalteparin: 1 dose/day
      ATC == "B01AB05" ~ doses, # Enoxaparin: 1 dose/day
      TRUE ~ NA_real_
    ),
    days_supply = ceiling(days_supply),
    supply_end = EDATUM + days_supply - 1
  )



DOAC %>%
  arrange(LopNr, EDATUM) %>%
  group_by(LopNr) %>%
  slice(1) %>%
  ungroup() %>%
  filter(ATC == "B01AF02", styrknum == 2.5)


# Calculate duration of each prescription
DOAC <- DOAC %>%
  mutate(
    tablets = ANTAL * antnum,
    days_supply = case_when(
      ATC == "B01AF01" ~ if_else(tablets <= 42, tablets / 2, 21 + (tablets - 42)), # Rivaroxaban: 2 tablets/day for 21 days, thereafter 1/day
      ATC == "B01AF02" & styrknum == 5 ~ if_else(tablets <= 28, tablets / 4, 7 + (tablets - 28) / 2), # Apixaban 5 mg: 4 tablets/day for 7 days, thereafter 2/day
      ATC == "B01AF02" & styrknum == 2.5 ~ tablets / 2, # Apixaban 2.5 mg: 2 tablets/day
      ATC %in% c("B01AF03", "B01AE07") ~ tablets / 2, # Edoxaban and Dabigatran: 2 tablets/day
      TRUE ~ NA_real_
    ),
    days_supply = ceiling(days_supply),
    supply_end = EDATUM + days_supply - 1
  )





get_treatment <- function(data) {

  data %>%
    filter(!is.na(EDATUM), !is.na(days_supply), days_supply > 0) %>%
    arrange(LopNr, EDATUM) %>%
    group_by(LopNr) %>%
    mutate(
      previous_end = as.Date(lag(cummax(as.numeric(supply_end))), origin = "1970-01-01"),
      new_episode = if_else(
        is.na(previous_end) | EDATUM > previous_end + 31,
        1L, 0L
      ),
      episode = cumsum(new_episode)
    ) %>%
    group_by(LopNr, episode) %>%
    mutate(
      treatment_start = min(EDATUM),
      treatment_end = max(supply_end)
    ) %>%
    ungroup()
}

DOAC <- get_treatment(DOAC)
library(dplyr)
library(arrow)

# Step 1: Combine all DOAC prescriptions
DOAC <- bind_rows(
  open_dataset(file.path(meds, "EVENTS_Rivaroxaban_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Apixaban_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Dabigatran_prescription.parquet")) %>% collect(),
  open_dataset(file.path(meds, "EVENTS_Edoxaban_prescription.parquet")) %>% collect()
) %>%
  filter(!(coalesce(ATC == "B01AF01" & styrknum == 2.5, FALSE))) %>%
  arrange(LopNr, EDATUM)

# Step 2: Calculate duration of each prescription
DOAC <- DOAC %>%
  mutate(
    tablets = ANTAL * antnum,
    days_supply = case_when(
      # Rivaroxaban: 2 tablets/day for 21 days, thereafter 1/day
      ATC == "B01AF01" ~
        if_else(tablets <= 42, tablets / 2, 21 + (tablets - 42)),

      # Apixaban: 4 tablets/day for 7 days, thereafter 2/day
      ATC == "B01AF02" ~
        if_else(tablets <= 28, tablets / 4, 7 + (tablets - 28) / 2),

      # Dabigatran and Edoxaban: 2 tablets/day
      ATC %in% c("B01AF03", "B01AF07") ~ tablets / 2,

      TRUE ~ NA_real_
    ),
    days_supply = ceiling(days_supply),
    supply_end = EDATUM + days_supply - 1
  )

# Step 3: Identify treatment episodes (30-day grace period)
get_treatment <- function(data) {

  data %>%
    filter(!is.na(EDATUM), !is.na(days_supply), days_supply > 0) %>%
    arrange(LopNr, EDATUM) %>%
    group_by(LopNr) %>%
    mutate(
      previous_end = lag(cummax(supply_end)),
      new_episode = if_else(
        is.na(previous_end) | EDATUM > previous_end + 31,
        1L, 0L
      ),
      episode = cumsum(new_episode)
    ) %>%
    group_by(LopNr, episode) %>%
    mutate(
      treatment_start = min(EDATUM),
      treatment_end = max(supply_end)
    ) %>%
    ungroup()
}

# Step 4: Add treatment episode information to each prescription
DOAC <- get_treatment(DOAC)

# Step 5: Inspect results
DOAC %>%
  select(LopNr, EDATUM, ATC, styrknum, days_supply, supply_end,
         episode, treatment_start, treatment_end) %>%
  arrange(LopNr, EDATUM)

















library(dplyr)

# DOAC: Calculate duration of each prescription
DOAC <- DOAC %>%
  mutate(
    tablets = Antalforpackning * packagesize,
    days_supply = case_when(
      # Rivaroxaban: 2 tablets/day for 21 days, thereafter 1/day
      Type == "B01AF01" ~
        if_else(tablets <= 42, tablets / 2, 21 + (tablets - 42)),

      # Apixaban: 4 tablets/day for 7 days, thereafter 2/day
      Type == "B01AF02" ~
        if_else(tablets <= 28, tablets / 4, 7 + (tablets - 28) / 2),

      # Dabigatran and Edoxaban: 2 tablets/day
      Type %in% c("B01AF03", "B01AF07") ~ tablets / 2,

      TRUE ~ NA_real_
    ),
    days_supply = ceiling(days_supply),
    supply_end = DOAC_date + days_supply - 1
  )

# LMWH: Calculate duration of each prescription
LMWH <- LMWH %>%
  mutate(
    doses = Antalforpackning * packagesize,
    days_supply = case_when(
      # Tinzaparin and Dalteparin: 1 dose/day
      Type %in% c("B01AB10", "B01AB04") ~ doses,

      TRUE ~ NA_real_
    ),
    days_supply = ceiling(days_supply),
    supply_end = LMWH_date + days_supply - 1
  )





library(dplyr)

# Step 1: Calculate duration of each prescription
TKI <- TKI %>%
  mutate(
    days_supply = ceiling(Antalforpackning * packagesize / daily_tablets),
    supply_end = TKI_date + days_supply - 1
  )

DOAC <- DOAC %>%
  mutate(
    days_supply = ceiling(Antalforpackning * packagesize / daily_tablets),
    supply_end = DOAC_date + days_supply - 1
  )

get_treatment <- function(data, date) {

  data %>%
    arrange(LopNr, {{ date }}) %>%
    group_by(LopNr) %>%
    mutate(
      previous_end = lag(cummax(supply_end)),
      new_episode = if_else(
        is.na(previous_end) | {{ date }} > previous_end + 31,
        1L, 0L
      ),
      episode = cumsum(new_episode)
    ) %>%
    group_by(LopNr, episode) %>%
    summarise(
      treatment_start = min({{ date }}),
      treatment_end = max(supply_end),
      .groups = "drop"
    )
}


library(dplyr)
library(lubridate)

# Step 1: Exclude patients aged <18 years
adsl <- VTE %>%
  mutate(excl = if_else(age < 18, 1L, 0L))

# Step 2: Combine anticoagulant prescriptions
AC <- bind_rows(
  DOAC %>% transmute(LopNr, AC_date = DOAC_date, treatment = "DOAC"),
  LMWH %>% transmute(LopNr, AC_date = LMWH_date, treatment = "LMWH")
)

# Step 3: Identify the earliest eligible baseline
baseline_dat <- adsl %>%
  filter(excl == 0) %>%
  select(LopNr, VTE_date) %>%
  inner_join(TKI, by = "LopNr", relationship = "many-to-many") %>%
  inner_join(AC, by = "LopNr", relationship = "many-to-many") %>%
  filter(
    AC_date >= TKI_date,
    AC_date <= TKI_date + 60,
    AC_date >= VTE_date,
    AC_date < as.Date("2025-10-01")
  ) %>%
  arrange(LopNr, AC_date, TKI_date) %>%
  group_by(LopNr) %>%
  slice(1) %>%
  ungroup() %>%
  rename(baseline = AC_date)

# Step 4: Exclude patients without eligible baseline
adsl <- adsl %>%
  left_join(
    baseline_dat %>% select(LopNr, baseline, treatment),
    by = "LopNr"
  ) %>%
  mutate(
    excl = case_when(
      excl == 0 & is.na(baseline) ~ 2L,
      TRUE ~ excl
    )
  )

# Step 5: Exclude previous DOAC or LMWH use
previous_AC <- adsl %>%
  filter(excl == 0) %>%
  select(LopNr, baseline) %>%
  inner_join(AC, by = "LopNr", relationship = "many-to-many") %>%
  filter(
    AC_date < baseline,
    AC_date >= baseline %m-% months(3)
  ) %>%
  distinct(LopNr) %>%
  mutate(previous_AC = 1L)

adsl <- adsl %>%
  left_join(previous_AC, by = "LopNr") %>%
  mutate(
    excl = case_when(
      excl == 0 & previous_AC == 1 ~ 3L,
      TRUE ~ excl
    )
  )

# Step 6: Exclude dialysis within 3 months prior to or at baseline
dialysis_excl <- adsl %>%
  filter(excl == 0) %>%
  select(LopNr, baseline) %>%
  inner_join(Dialysis, by = "LopNr", relationship = "many-to-many") %>%
  filter(
    dialysis_date >= baseline %m-% months(3),
    dialysis_date <= baseline
  ) %>%
  distinct(LopNr) %>%
  mutate(dialysis = 1L)

adsl <- adsl %>%
  left_join(dialysis_excl, by = "LopNr") %>%
  mutate(
    excl = case_when(
      excl == 0 & dialysis == 1 ~ 4L,
      TRUE ~ excl
    )
  )

# Step 7: Exclude strong CYP3A4 inhibitor use
cyp_excl <- adsl %>%
  filter(excl == 0) %>%
  select(LopNr, baseline) %>%
  inner_join(CYP3A4, by = "LopNr", relationship = "many-to-many") %>%
  filter(
    CYP3A4_date >= baseline %m-% months(3),
    CYP3A4_date <= baseline
  ) %>%
  distinct(LopNr) %>%
  mutate(CYP3A4 = 1L)

adsl <- adsl %>%
  left_join(cyp_excl, by = "LopNr") %>%
  mutate(
    excl = case_when(
      excl == 0 & CYP3A4 == 1 ~ 5L,
      TRUE ~ excl
    )
  )

# Step 8: Summary of sequential exclusions
adsl %>%
  count(excl)

# Step 9: Final analysis population
adsl_FAS <- adsl %>%
  filter(excl == 0)

### Read in Total population register
SCB <- read_parquet("/safe/data/SCARP/DATA_2024/PARQUET/ar_lev_grunduppgifter.parquet") %>% rename(LopNr=lopnr)

### Population
Pop <- open_dataset(file.path(der,"FirstTimeOutcome_HF.parquet")) %>% collect() %>%
        transmute(LopNr=LopNr, baseline = INDATUM, lisayr = year(baseline)) %>%
        left_join(SwedeHF, by = "LopNr") %>%
        mutate(in_SwedeHF = as.integer(LopNr %in% SwedeHF$LopNr))


### Exclusion criteria:
Pop1 <- Pop %>% mutate(excl = if_else(is.na(baseline),1,0))

Pop2 <- Pop1 %>% mutate(excl = if_else(excl == 0 & (baseline<as.Date("2009-01-01") | baseline==as.Date("2023-12-31")),2,excl))

Pop3 <- Pop2 %>% left_join(open_dataset(file.path(der,"all_cause_death.parquet"))%>% collect(),by="LopNr") %>%
        mutate(deathdt = case_when(
            deathdt <= baseline & !is.na(befdoddtm) ~ befdoddtm, # These are death dates with wrong format which is why heart failure can happen afterwards
            deathdt <= baseline & is.na(befdoddtm) ~ baseline,  # We assume if it is available in Riksvikt or assume it is the same day as baseline
            TRUE ~ deathdt),
        excl=case_when(excl==0 & (deathdt == baseline)~3,TRUE~excl))

excl_aterpnr <- SCB %>% filter(aterpnr == 1) %>% distinct(LopNr) %>% mutate(excl_aterpnr = 2) %>% select(LopNr,excl_aterpnr)

Pop4 <- Pop3 %>% left_join(excl_aterpnr,by ="LopNr") %>%
        ungroup() %>%
        mutate(excl = case_when(excl == 0 & excl_aterpnr ==4 ~ excl_aterpnr,TRUE~excl)) %>%
        select(-excl_aterpnr)

save(Pop4, file = file.path(der,"Pop4.RData"))
load(file.path(der,"Pop4.RData"))

Pop5 <- Pop4 %>%
        left_join(open_dataset(file.path(der,"CIED_exclusion.parquet")) %>% collect(),by="LopNr") %>%
        mutate(prior_cied=ifelse(!is.na(EVENTDATE_CIED) & EVENTDATE_CIED <=baseline, 1, 0),
               excl=case_when(excl==0 & prior_cied==1~5,TRUE~excl))

Pop6 <- Pop5 %>%
        mutate(excl=case_when(excl==0 & !REASON_CIED %in% c("Primärimplantation","N/A",NA)==1~6,TRUE~excl),
               FAS = if_else(excl == 0,1,0)) %>%
        select(-c(prior_cied,EVENTDATE_CIED,REASON_CIED,befdoddtm))


##save to continue later
save(Pop6, file = file.path(der,"Pop6.RData"))
load(file.path(der,"Pop6.RData"))









get_disease(Disease = "ALL", ICD10 = "^C91") # Acute lymphatic leukemia
get_disease(Disease = "Breast_cancer", ICD10 = "^C50") # Breast cancer
get_disease(Disease = "CML", ICD10 = "^C921") # Chronic myeloid leukemia
get_disease(Disease = "GIST", ICD10 = "^C1[6-8]") # Gastrointestinal stromal tumor
get_disease(Disease = "Kidney_cancer", ICD10 = "^C64") # Kidney cancer
get_disease(Disease = "Lung_cancer", ICD10 = "^C34") # Lung cancer
get_disease(Disease = "pNET", ICD10 = "^C254") # Pancreatic neuroendocrine tumor
get_disease(Disease = "Thyroid_cancer", ICD10 = "^C73") # Thyroid cancer
get_disease(Disease = "Urothelial_cancer", ICD10 = "^C659|^C669|^C67|^C680|^C688|^C689") # Urothelial cancer





# Diabetes
get_disease(Disease = "Diabetes", ICD10 = "^E1[0-4]")

# Renal disease
get_disease(Disease = "Renal_disease", ICD10 = "^N1[7-9]")

# Liver disease
get_disease(Disease = "Liver_disease", ICD10 = "^K7[0-7]")

# Ischemic heart disease
get_disease(Disease = "Ischemic_heart_disease", ICD10 = "^I2[0-5]")

# Peripheral arterial disease
get_disease(Disease = "Peripheral_arterial_disease", ICD10 = "^I739|^I70")

# Ischemic stroke
get_disease(Disease = "Ischemic_stroke", ICD10 = "^I63|^I64|^G459")

# Hemorrhagic stroke
get_disease(Disease = "Hemorrhagic_stroke", ICD10 = "^I6[0-2]")

# Systemic connective tissue disorders
get_disease(Disease = "Systemic_connective_tissue_disorders", ICD10 = "^M3[0-6]")

# Inflammatory bowel disease
get_disease(Disease = "Inflammatory_bowel_disease", ICD10 = "^K5[0-2]")















# Intracranial bleeding, not traumatic
get_disease(Disease = "Nontraumatic_intracranial_bleeding", ICD10 = "^I60|^I61|^I62")

# Intracranial bleeding, traumatic
get_disease(Disease = "Traumatic_intracranial_bleeding", ICD10 = "^S064|^S065|^S066")

# Hemopericardium
get_disease(Disease = "Hemopericardium", ICD10 = "^I312")

# Gastrointestinal bleeding
get_disease(Disease = "Gastrointestinal_bleeding", ICD10 = "^I850|^I983|^K250|^K252|^K254|^K256|^K260|^K262|^K264|^K266|^K270|^K272|^K274|^K276|^K280|^K282|^K284|^K286|^K290|^K625|^K922|^K920|^K921|^K661")

# Bleeding in the airways
get_disease(Disease = "Airway_bleeding", ICD10 = "^R04|^R58|^J942|^S271|^H922")

# Anemia due to bleeding
get_disease(Disease = "Anemia_due_to_bleeding", ICD10 = "^D629")

# Genital bleeding
get_disease(Disease = "Genital_bleeding", ICD10 = "^N920|^N921|^N922|^N923|^N924|^N93|^N950|^N421|^N501A")

# Hematuria
get_disease(Disease = "Hematuria", ICD10 = "^R31")

# Bleeding in the eye
get_disease(Disease = "Eye_bleeding", ICD10 = "^H113|^H313|^H356|^H431|^H450")







# LISA
path_0 <- "/safe/data/Original data 2026/SCB/Leverans"
path_1 <- "/safe/data/Original data 2026/PARQUET"

for (i in 2009:2022) {

  lisa <- read_sas(file.path(path_0, paste0("ar_lev_lisa_", i, ".sas7bdat")))

  write_parquet(lisa, file.path(path_1, paste0("ar_lev_lisa_", i, ".parquet")))

  rm(lisa)
  gc()

}




if (!is.null(ICD_exclude)) {
  if (position != 1) {
    paste0("filter(!if_any(.cols=c(HDIA,DIA1:DIA", number_of_columns,
           "),.fns=list(G=~grepl('", ICD_exclude, "',.x)))) %>%")
  } else {
    paste0("filter(!grepl('", ICD_exclude, "',HDIA)) %>%")
  }
} else {
  ""
},



death_compare <- death %>%
  select(LopNr, DODSDAT, deathdt) %>%
  rename(
    DODSDAT_death = DODSDAT,
    deathdt_death = deathdt
  ) %>%
  left_join(
    death_avi %>%
      select(LopNr, DODSDAT, deathdt) %>%
      rename(
        DODSDAT_avi = DODSDAT,
        deathdt_avi = deathdt
      ),
    by = "LopNr"
  )


death_compare %>%
  filter(!is.na(deathdt_avi),
         deathdt_death != deathdt_avi)




library(dplyr)
library(lubridate)

# Dates affected by format_DODSDAT
death_changed <- death %>%
  filter(deathdt != ymd(DODSDAT, quiet = TRUE) |
         (is.na(ymd(DODSDAT, quiet = TRUE)) & !is.na(deathdt)))

avi_changed <- death_avi %>%
  filter(deathdt != ymd(DODSDAT, quiet = TRUE) |
         (is.na(ymd(DODSDAT, quiet = TRUE)) & !is.na(deathdt)))

# Patients with modified dates in the main register but not AVI
death_changed_only <- death_changed %>%
  anti_join(avi_changed, by = "LopNr")

nrow(death_changed_only)




library(arrow)
library(dplyr)
library(lubridate)

# Convert DODSDAT to Date, replacing unknown month/day with 01.
format_DODSDAT <- function(x) {
  x <- as.character(x)
  x <- trimws(x)

  date_string <- case_when(
    nchar(x) == 8 & substr(x, 5, 8) == "0000" ~ paste0(substr(x, 1, 4), "0101"),
    nchar(x) == 8 & substr(x, 7, 8) == "00" ~ paste0(substr(x, 1, 6), "01"),
    nchar(x) == 6 & substr(x, 1, 2) %in% c("19", "20") ~ paste0(x, "01"),
    TRUE ~ x
  )

  ymd(date_string, quiet = TRUE)
}

# Main death register
death <- open_dataset("/safe/data/SCARP/DATA_2024/PARQUET/ut_r_dors_19328_2023.parquet") %>%
  collect()

nrow(death)
n_distinct(death$LopNr)

d1 <- death %>%
  mutate(deathdt = format_DODSDAT(DODSDAT))

# Supplementary death register (AVI)
death_avi <- open_dataset("/safe/data/SCARP/DATA_2024/PARQUET/ut_r_dors_avi_19328_2023.parquet") %>%
  collect()

da1 <- death_avi %>%
  mutate(deathdt = format_DODSDAT(DODSDAT))

# Keep all main-register records; add AVI records only for new patients.
all_cause_death <- bind_rows(
  d1 %>% select(LopNr, DODSDAT, deathdt),
  da1 %>% anti_join(d1, by = "LopNr") %>% select(LopNr, DODSDAT, deathdt)
)

# Retain first recorded death on or after 2000-01-01 for each patient.
all_cause_death <- all_cause_death %>%
  filter(!is.na(deathdt), deathdt >= as.Date("2000-01-01")) %>%
  arrange(LopNr, deathdt) %>%
  group_by(LopNr) %>%
  slice(1) %>%
  ungroup()

n_distinct(all_cause_death$LopNr)
nrow(all_cause_death)

write_parquet(
  all_cause_death,
  "/safe/data/APNC/SCARP-HFDEVICE/Data/Derived/all_cause_death.parquet"
)












" mutate(across(c(hdia,DIA1:DIA",number_of_columns,"),~sub('\\\\\\.|\\\\\\,|\\\\\\-','','.',.)) %>%
        filter(AR>=",Year_min,",AR<",Year_max,")  %>%",
        if_else(position!=1,
                paste0("filter(if_any(.cols=c(hdia,DIA1:DIA",number_of_columns,"),.fns=list(G=~grepl('",ICD_string,"',.x)))) %>%"),
                paste0("filter(grepl('",ICD_string,"',hdia)) %>%")),
        "group_by(LopNr) %>%",
        dplyr::case_when(REG=="SV"~"arrange(UTDATUM,.by_group = TRUE) %>%",
                         REG=="OV"~"mutate(UTDATUM=INDATUM) %>% arrange(INDATUM,.by_group = TRUE) %>%",
                         REG=="CDR"~"mutate(INDATUM=UTDATUM) %>% arrange(UTDATUM,.by_group = TRUE) %>%",
                         REG=="VPR"~"mutate(INDATUM=UTDATUM) %>% arrange(UTDATUM,.by_group = TRUE) %>%",
                         TRUE~"arrange(INDATUM,.by_group = TRUE) %>%"),
        "collect() %>%",# mutate(id=row_number()) %>%
        #filter(id==1) %>%",
        dplyr::case_when(REG=="SV"~"select(LopNr,AR,INDATUM,UTDATUM,hdia,DIA1:DIA30,ALDER,KON) %>%",
                         REG=="OV"~"select(LopNr,AR,INDATUM,UTDATUM,hdia,DIA1:DIA30) %>%",
                         REG=="CDR"~"select(LopNr,AR,UTDATUM,INDATUM,hdia,DIA1:DIA48) %>%",
                         REG=="VPR"~"select(LopNr,AR,UTDATUM,INDATUM,hdia,DIA1:DIA6) %>%",
                         TRUE~"select(LopNr,AR,INDATUM,UTDATUM,hdia,DIA1:DIA30,ALDER,KON) %>%"),
        "identity()")
apa <- eval(parse(text=apa_string)) %>% mutate(Diagnose=Name,Reg=REG)

return(apa)
}








Create_dataset <- function(Disease = "Overall_CV_Diag",
                           Register = c("CDR", "SV", "OV"),
                           Position = 2,
                           YearICD10 = c(1996, Inf),
                           ICD10 = "^I|^Q2[2-3]|^Q243<^Q244") {

  X_list <- list()

  # INPATIENT ICD10
  if ("SV" %in% Register) {
    X_list[["SV"]] <- parquet_get_dia(
      "SV", Name = Disease, ICD_string = ICD10,
      position = Position,
      Year_min = YearICD10[1], Year_max = YearICD10[2]
    ) %>%
      select(-ALDER, -KON)
  }

  # OUTPATIENT ICD10
  if ("OV" %in% Register) {
    X_list[["OV"]] <- parquet_get_dia(
      "OV", Name = Disease, ICD_string = ICD10,
      position = Position,
      Year_min = YearICD10[1], Year_max = YearICD10[2]
    ) %>%
      select(-ALDER, -KON)
  }

  # CDR ICD10
  if ("CDR" %in% Register) {
    X_list[["CDR"]] <- parquet_get_dia(
      "CDR", Name = Disease, ICD_string = ICD10,
      position = Position,
      Year_min = YearICD10[1], Year_max = YearICD10[2]
    ) %>%
      select(-ALDER, -KON)
  }

  # Combine registers and keep first diagnosis
  X <- bind_rows(X_list) %>%
    group_by(LopNr) %>%
    arrange(INDATUM, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(LopNr, INDATUM, Diagnose)

  # Add disease-specific variables
  X[, paste0("INDATUM_", Disease)] <- X$INDATUM
  X[, Disease] <- 1

  # Save only final dataset
  X %>%
    write_parquet(
      paste0(
        "/safe/data/APNC/SCARLET morbidity/Data/Derived/FirstTimeOutcome_",
        Disease, ".parquet"
      )
    )

  rm(X_list, X)
  gc()
}







library(reticulate)

py_install("pyreadstat")

pyreadstat <- import("pyreadstat")

result <- pyreadstat$read_sas7bdat(
  sas_path,
  metadataonly = TRUE
)

result[[2]]$number_rows




library(sas7bdat)

sas_path <- file.path(
  path_0,
  "ut_r_par_ov_ut1_34759_2025.sas7bdat"
)

header <- read.sas7bdat.header(sas_path)

header



library(arrow)
library(dplyr)

parquet_path <- file.path(path_1, "ut_r_par_ov_ut1_34759_2025")

ds <- open_dataset(parquet_path, format = "parquet")

ds %>%
  summarise(n = n()) %>%
  collect()


proc sql;
    select count(*) as N
    from original.ut_r_par_ov_ut1_34759_2025;
quit;



data data;
    set adsl(obs=10000);
    where FAS=1;

    if CIEDcomprisk = "Censored" then CIEDcomprisk_n = 0;
    else if CIEDcomprisk = "Event" then CIEDcomprisk_n = 1;
    else if CIEDcomprisk = "Death" then CIEDcomprisk_n = 2;

    drop CIEDcomprisk;
    rename CIEDcomprisk_n = CIEDcomprisk;

    keep LopNr excl CIEDcomprisk_n CIEDcompriskyrs
         EVENTDATE_CIED deathdt baseline cens;
run;
library(parquetize)
library(arrow)
library(dplyr)

path_0 <- "/safe/data/Original data 2026/Komplettering öppenvård"
path_1 <- "/safe/data/Original data 2026/PARQUET"

max_mem <- 2000

files <- c(
  "ut_r_par_ov_ut1_34759_2025",
  "ut_r_par_ov_ut2_34759_2025"
)

# Step 1: Convert SAS files to Parquet
for (file in files) {

  table_to_parquet(
    path_to_file = file.path(path_0, paste0(file, ".sas7bdat")),
    path_to_parquet = file.path(path_1, paste0(file, ".parquet")),
    max_memory = max_mem
  )

}

# Step 2: Combine Parquet files
dat <- bind_rows(
  open_dataset(file.path(path_1, paste0(files[1], ".parquet"))) %>% collect(),
  open_dataset(file.path(path_1, paste0(files[2], ".parquet"))) %>% collect()
)

# Step 3: Save combined Parquet file
write_parquet(
  dat,
  file.path(path_1, "ut_r_par_ov_34759_2025.parquet")
)



# Create time-dependent dataset for heart transplantation
adtte_trans <- tmerge(
  data1 = adsl,
  data2 = adsl,
  id = LopNr,
  tstop = censyrs,
  trans_event = event(TRANSYRS, TRANS)
)

adtte_trans <- tmerge(
  data1 = adtte_trans,
  data2 = adsl,
  id = LopNr,
  PM_td = tdc(PMyrs, init = 0),
  ICD_td = tdc(ICDyrs, init = 0),
  CRT_td = tdc(CRTyrs, init = 0),
  CIED_td = tdc(CIEDyrs, init = 0)
) %>%
  mutate(person_time = tstop - tstart)





line "PY = Person-years; IR = Incidence rate; HR = Hazard ratio.";
line "Incidence rates are expressed per 100 person-years. Deaths are attributed to the exposure status at the time of death. Patients may contribute person-time to both exposure groups.";
line "Adjusted hazard ratios are adjusted for sex and age.";


	value visitftwo 1 = "Visit 1 (Day 0)"
	             2 = "¤{style[font_style=italic just=left]Secondary Endpoint 2}¤nVisit 2 (1 month ±14 days)"
				 3 = "Visit 3 (3 months ±14 days)"
				 4 = "¤{style[font_style=italic just=left]Secondary Endpoint 2}¤nVisit 4 (6 months ±14 days)";



# Unadjusted time-dependent Cox regression
fit <- coxph(
  Surv(tstart, tstop, death_event) ~ PM_td,
  data = adtte
)

summary(fit)

# HR, 95% CI and p-value
broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
fit_adj <- coxph(
  Surv(tstart, tstop, death_event) ~ PM_td + sex + age + Birth_region,
  data = adtte
)

summary(fit_adj)

broom::tidy(fit_adj, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "PM_td")


tdc_label <- attr(data[[tdc]], "label")
if (is.null(tdc_label)) tdc_label <- tdc

result <- data.frame(
  outcome = outcome,
  tdc = tdc_label,
  group = c("No", "Yes"),
  events = desc$events[idx],
  PY = round(desc$PY[idx], 1),
  IR = round(desc$IR[idx], 2),
  HR = c("Reference", HR),
  p_value = c("", p)
)









# All-cause mortality
result_death <- purrr::map(
  c("PM_td", "ICD_td", "CRT_td", "CIED_td"),
  ~ td_cox_table(
    data = adtte,
    tdc = .x,
    outcome = "death_event",
    adj = c("sex", "age", "Birth_region")
  )
) %>%
  bind_rows()


# Cardiovascular mortality
result_cv <- purrr::map(
  c("PM_td", "ICD_td", "CRT_td", "CIED_td"),
  ~ td_cox_table(
    data = adtte,
    tdc = .x,
    outcome = "cv_death_event",
    adj = c("sex", "age", "Birth_region")
  )
) %>%
  bind_rows()


# Sudden cardiac death
result_scd <- purrr::map(
  c("PM_td", "ICD_td", "CRT_td", "CIED_td"),
  ~ td_cox_table(
    data = adtte,
    tdc = .x,
    outcome = "scd_death_event",
    adj = c("sex", "age", "Birth_region")
  )
) %>%
  bind_rows()


# Combine all results
result <- bind_rows(result_death, result_cv, result_scd)






outcome_label <- attr(data[[outcome]], "label")
if (is.null(outcome_label)) outcome_label <- outcome


td_dat <- tmerge(
  data1 = adsl,
  data2 = adsl,
  id = LopNr,
  tstop = censyrs,
  death_event = event(deathyrs, death),
  cv_death_event = event(deathyrs, CVdeath == 1),
  scd_death_event = event(deathyrs, SCDdeath == 1)
)



result <- purrr::map(
  c("PM_td", "ICD_td", "CRT_td", "CIED_td"),
  ~ td_cox_table(
      data = td_dat,
      tdc = .x,
      outcome = "death_event",
      adj = c("sex", "Birth_region")
    )
) %>%
  bind_rows()



# Formats
pval_func <- function(pval) {
  case_when(
    is.na(pval)    ~ "",
    pval > 0.1     ~ sprintf("%.2f", pval),
    pval >= 0.0001 ~ sprintf("%.3f", pval),
    TRUE           ~ "<.0001"
  )
}

# Convert Yes/No to 1/0
yes_no <- function(x) {
  case_when(
    x == "YES" ~ 1,
    x == "NO"  ~ 0,
    TRUE       ~ NA_real_
  )
}

ynf <- function(x) {
  factor(x, levels = 0:1, labels = c("No", "Yes"))
}

prep_func_new <- function(var) {
  var <- na_if(var, ".A")
  var <- na_if(var, ".N")
  var <- na_if(var, ".E")
  var <- na_if(var, "")
}

# Bounded complementary log-log link
bcloglog <- function(edge = 0.5) {
  link <- make.link("cloglog")
  link$linkfun <- function(mu) {
    x <- pmax(edge, pmin(mu, 1 - edge))
    log(-log(1 - x))
  }
  link$name <- "bcloglog"
  link
}

# Format estimates and 95% CI
fmt_ci <- function(est, se) {
  sprintf("%.2f (%.2f - %.2f)",
          100 * est,
          100 * (est - 1.96 * se),
          100 * (est + 1.96 * se))
}

# First age where CIF exceeds threshold
first_age_above <- function(df, threshold) {
  ages <- df$time[df$cif > threshold]
  if (length(ages) == 0) NA_real_ else min(ages)
}

cif_at_time <- function(df, time_point) {
  idx <- which(df$time >= time_point)
  if (length(idx) == 0) NA_real_ else df$cif[idx[1]]
}


td_cox_table <- function(
    data,
    tdc,
    outcome,
    adj = NULL,
    tstart = "tstart",
    tstop = "tstop"
) {
  
  # Descriptive statistics
  desc <- data %>%
    group_by(.data[[tdc]]) %>%
    summarise(
      events = sum(.data[[outcome]], na.rm = TRUE),
      PY = sum(person_time),
      IR = 100 * events / PY,
      .groups = "drop"
    )
  
  # Helper function for Cox regression
  get_cox <- function(vars) {
    
    formula <- reformulate(
      vars,
      response = paste0(
        "Surv(", tstart, ", ", tstop, ", ", outcome, ")"
      )
    )
    
    coxph(formula, data = data) %>%
      broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == tdc)
  }
  
  # Unadjusted Cox regression
  s <- get_cox(tdc)
  
  HR <- fmt_ci(s$estimate, s$conf.low, s$conf.high)
  p <- pval_func(s$p.value)
  
  # Combine results
  idx <- match(0:1, desc[[tdc]])
  
  result <- data.frame(
    outcome = outcome,
    group = paste(tdc, c("No", "Yes")),
    events = desc$events[idx],
    PY = round(desc$PY[idx], 1),
    IR = round(desc$IR[idx], 2),
    HR = c("Reference", HR),
    p_value = c("", p)
  )
  
  # Adjusted Cox regression
  if (!is.null(adj)) {
    
    s_adj <- get_cox(c(tdc, adj))
    
    result$HR_adj <- c(
      "Reference",
      fmt_ci(s_adj$estimate, s_adj$conf.low, s_adj$conf.high)
    )
    
    result$p_value_adj <- c("", pval_func(s_adj$p.value))
  }
  
  return(result)
}








td_cox_table(
  data = td_dat,
  tdc = "TD1",
  outcome = "death",
  adj = c("age", "sex")
)


set.seed(123)

# Generate patient-level data
adsl <- data.frame(
  id = 1:100,
  age = round(rnorm(100, 65, 10)),
  sex = sample(0:1, 100, replace = TRUE),
  deathyrs = runif(100, 1, 10),
  death = rbinom(100, 1, 0.4),
  TD1yrs = runif(100, 0, 10),
  TD2yrs = runif(100, 0, 10)
)

# Create time-dependent dataset
td_dat <- tmerge(
  data1 = adsl, data2 = adsl, id = id,
  tstop = deathyrs,
  death = event(deathyrs, death)
)

td_dat <- tmerge(
  data1 = td_dat, data2 = adsl, id = id,
  TD1 = tdc(TD1yrs, init = 0),
  TD2 = tdc(TD2yrs, init = 0)
) %>%
  select(id, tstart, tstop, death, TD1, TD2, age, sex) %>%
  mutate(person_time = tstop - tstart)

# Run time-dependent Cox regression
result <- bind_rows(
  td_cox_table(td_dat, "TD1", "death", adj = c("age", "sex")),
  td_cox_table(td_dat, "TD2", "death", adj = c("age", "sex"))
)

print(result)

fmt_ci <- function(est, lower, upper) {
  sprintf("%.2f (%.2f - %.2f)", est, lower, upper)
}




### Create empty data frames to fill with values
Pop <- open_dataset("/safe/data/APNC/SCARP-HFDEVICE/Data/Derived/FirstTimeOutcome_HF.parquet") %>% collect() %>% select(LopNr)
lisa_pop <- data.frame()

for (i in 2009:2022) {
  lisa<- open_dataset(paste0("/safe/data/SCARP/DATA_2024/PARQUET/ar_lev_lisa_", i, ".parquet")) %>% collect()
  colnames(lisa)<- tolower(colnames(lisa))
  if (i>=2005 & i<2019){
    lisa<- lisa %>% select(lopnr, dispinkfam04, dispink04, sun2000niva_old) %>% rename("income_ind"="dispink04" , "education_level"="sun2000niva_old") %>% rename(LopNr=lopnr)
  } else if (i >=2019 ){
    lisa<- lisa %>% select(lopnr, dispinkfam04, dispink04, sun2020niva_old) %>% rename("income_ind"="dispink04" , "education_level"="sun2020niva_old") %>% rename(LopNr=lopnr)
  }

  lisa_i<- Pop %>% left_join(lisa, by ="LopNr") %>% mutate(lisayr=i)
  lisa_pop <- rbind(lisa_pop, lisa_i)

  rm(lisa_i, lisa)
  print(paste("year", i, "is done"))# print to see progress
}

# Make education level numeric and remove the stars
lisa_pop <- lisa_pop %>% mutate(education_level = as.numeric(na_if(education_level, "*")))

lisa_pop_LOCF<- lisa_pop %>%
  group_by(LopNr) %>%
  fill(education_level, income_ind, .direction="down") %>%
  ungroup() %>%
  arrange(LopNr, lisayr) %>%
  rename("income"="income_ind", "education_level" = "education_level") %>%
  select(LopNr,income,education_level,lisayr)

KPI_2009_2022 <- data.frame(lisayr = 2009:2022,
                           KPI = c(299.66, 303.46, 311.43, 314.20, 314.06, 313.49, 313.35, 316.43, 322.11, 328.40, 334.26, 335.92, 343.19, 371.91)) %>%
  mutate(KPI_factor = 299.66 / KPI)

lisa_pop_LOCF <- lisa_pop_LOCF %>%
  left_join(KPI_2009_2022, by="lisayr") %>%
  mutate(income_KPI_adj=round(income*KPI_factor))

lisa_pop_LOCF %>% write_parquet("/safe/data/APNC/SCARP-HFDEVICE/Data/Derived/LISA.parquet")







######################################
#####  Specific CVDs             #######
######################################
# CV="CV"
# ICD_string="^A28|^I50|^I110|^I130|^I132"

spec_cv_fun<- function(CV="CV", ICD_string="^XXX"){

  X_SV<- parquet_get_dia("SV",Name=CV, ICD_string = ICD_string , position = 2,Year_min = 1970,Year_max = 2025)
  X_OV<- parquet_get_dia("OV",Name=CV, ICD_string = ICD_string, position = 2,Year_min = 1987,Year_max = 2025)

  X_SV<- X_SV %>% select(-ALDER,-KON)
#  colnames(X_SV)==colnames(X_OV) #kontroll
  X<- rbind(X_SV, X_OV)
  #Only need first occurence
  X<- X %>% group_by(LopNr) %>% arrange(INDATUM) %>% slice(1) %>% select(LopNr, INDATUM, Diagnose)

  X[, paste0("INDATUM_", CV)]<- X$INDATUM # keep both with same name for all data sets and a specific date name for now- not sure how to
  X[, paste0(CV)]<- 1 ## All have had an event, filtered data
  #excl_downs_syndome<- CCS_kohort %>% inner_join(y = excl_downs_syndome, by = "LopNr") #161 ## TA BORT DETTA! BEHÖVS ÄVEN FÖR KONTROLLER

  X %>% write_parquet(paste0("/safe/data/APNC/SCARLET morbidity/Data/Derived/FirstTimeOutcome_",CV,".parquet"))

  rm(X_OV, X_SV, X)

}

#Testa
spec_cv_fun(CV="HF", ICD_string = "^A28|^I50|^I110|^I130|^I132") #heart failure
spec_cv_fun(CV="MI", ICD_string = "^410|^I21")# myocardial infarction
spec_cv_fun(CV="angina", ICD_string = "^411B|^413|^I20") #angina pectoris
spec_cv_fun(CV="valvular", ICD_string = "^394|^395|^396|^397|^424|^I05|^I06|^I07|^I08|^I091|^I34|^I35|^I36|^I37|^I38") #valvular disease
spec_cv_fun(CV="card_arrest", ICD_string = "^427E|^I427F|^I46") #cardiac arrest
spec_cv_fun(CV="AFib", ICD_string = "^427D|^I48") #atrial fibrillation of flutter
## INSERT CONDUCTION DISEASE WHEN DONE!
spec_cv_fun(CV="stroke", ICD_string = "^430|^431|^434|^436|^I60|^I61|^I63") #stroke
spec_cv_fun(CV="PAD", ICD_string = "^440C|^444C|^I702|^I739") #pad
spec_cv_fun(CV="aortic", ICD_string = "^441A|^441B|^441D|^I710|^I711|^I713|^I715") #aortic disease
spec_cv_fun(CV="VTE", ICD_string = "^415B|^416W|^451|^I26|^I80") #venous thromboembolism

spec_cv_fun(CV="CM", ICD_string = "^425|^I42|^I43") #cardiomyopathy
spec_cv_fun(CV="HCM", ICD_string = "^425B|^I421|^I422") #cardiomyopathy: hypertropic
spec_cv_fun(CV="DCM", ICD_string = "^425E|^425F|^I420|^I426") #cardiomyopathy: dilated
spec_cv_fun(CV="RCM", ICD_string = "^I425") #cardiomyopathy: restrictive
spec_cv_fun(CV="OCM", ICD_string = "^425A|^425D|^425E|^425H|^425X|^425W|^I423|^I424|^I427|^I428|^I429|^I430|^I431|^I432|^I438|^O903") #cardiomyopathy: other
spec_cv_fun(CV="FCM", ICD_string = "^I427") #cardiomyopathy: by farmacological agent

spec_cv_fun(CV="IHD", ICD_string = "^391A|^393|^420|^423|^I30|^I31|^I32|^I010|^391C|^422|^I40|^I41|^I012|^391B|^421|^I011|^I33|^I398") #inflammatory heart disease
spec_cv_fun(CV="IHD_pericard", ICD_string = "^391A|^393|^420|^423|^I30|^I31|^I32|^I010") #pericarditis
spec_cv_fun(CV="IHD_myocard", ICD_string = "^391C|^422|^I40|^I41|^I012") #myocarditis
spec_cv_fun(CV="IHD_endocard", ICD_string = "^391B|^421|^I011|^I33|^I398") #endocarditis

spec_cv_fun(CV="Hypertension", ICD_string="^402|^I10|^I11|^I12|^I13|^I14|^I15") #Hypertension























library(haven)
library(arrow)
library(dplyr)

path_0 <- "/safe/data/Original data 2026/PARQUET"

dat <- bind_rows(
  read_sas(file.path(path_0, "ut_r_par_ov_ut1_34759_2025.sas7bdat")),
  read_sas(file.path(path_0, "ut_r_par_ov_ut2_34759_2025.sas7bdat"))
)

write_parquet(dat, file.path(path_0, "ut_r_par_ov_34759_2025.parquet"))




library(arrow)
library(dplyr)

path <- "your/path"

dat <- bind_rows(
  open_dataset(file.path(path, "ut_uttag1_r_par_sv_34759_2025.parquet")) %>% collect(),
  open_dataset(file.path(path, "ut_uttag2_r_par_sv_34759_2025.parquet")) %>% collect()
)

write_parquet(dat, file.path(path, "ut_combined.parquet"))




proc sql;
    create table data as
    select a.*,
           case when a.match_id in (
               select distinct match_id
               from der.adsl_1
               where LopNr in (select LopNr from B)
                 and not missing(match_id)
           ) then 1 else 0 end as B
    from der.adsl_1 as a;
quit;
df %>%
  filter(if_any(c(hdia, DIA1, DIA2, DIA3, DIA4),
                ~ substr(., 4, 4) == "A"))
