adsl <- adsl %>%
  mutate(
    deathdt = if_else(
      deathdt <= baseline & cause_of_death == "CV death",
      baseline,
      deathdt
    )
  )




*** get N as a macro variable;
proc sql;
select Placebo, Dexamethasone
into :N_Placebo, :N_Dexamethasone
from T1 where A='';
quit;

%put &N_Placebo. &N_Dexamethasone.;

data T2;
	set T2(rename=(Placebo=Placebo_tmp
	               Dexamethasone=Dexamethasone_tmp));
	if A ='' then delete;

	Placebo = tranwrd(Placebo_tmp,";","¤n");
	Dexamethasone = tranwrd(Dexamethasone_tmp,";","¤n");

	drop Placebo_tmp Dexamethasone_tmp Overall_tmp;
run;

data hj;
	length A $200;
	s = 0.5;	A = "¤S={just = left font_weight = bold font_style = italic}Corticosteroid before study start";	output;
	s = 9.5;	A = "¤S={just = left font_weight = bold font_style = italic}Corticosteroids during study intervention";	output;
run;

data T2;
    set hj T2;
    if substr(A,1,1) = ' ' then A = cats("¤S={leftmargin=10}", strip(A));
	format A;
run;

proc sort data=T2;
    by s;
run;

options nodate nonumber leftmargin=2.5cm rightmargin=2.5cm topmargin=2.5cm bottommargin=2.5cm;
ods escapechar ="¤";
ods rtf file="&projroot.\Output\TableHeader_02_&rdat8..rtf" style=styles.statistical bodytitle;
title;
ods rtf text= "Table 2. Corticosteroid exposure (PP population)";
ods rtf close;


ods escapechar ="¤";
ods rtf file="&projroot.\Output\Table_02_&rdat8..rtf" style=styles.statistical bodytitle;
proc report data= T2 nowd split="#" style(header)={background=CXECEDEC};
    column A Placebo Dexamethasone D;
    define A/ "Variable" style= {just=left cellwidth=4.5cm font_size=1};
    define Placebo/ "Placebo¤n &N_Placebo." style= {just=center cellwidth=3cm font_size=1};
    define Dexamethasone/ "Dexamethasone¤n &N_Dexamethasone." style= {just=center cellwidth=3cm font_size=1};
    define D/ "p-value" style= {just=center cellwidth=1.2cm font_size=1};

    compute after / style={just=left font_size=1 backgroundcolor=CXECEDEC };
    line "Categorical variables are presented as count(%) and continuous variables are presented by median [range] and number of non-missing.";
	line "Categorical variables were compared using Fisher’s exact test, and continuous variables were compared using the Mann–Whitney U test.";
    line "SAS version 9.4";
    endcomp;
run;
ods rtf close;














attr(adsl$bl_MI, "label")           <- "Myocardial infarction"
attr(adsl$bl_angina, "label")       <- "Angina pectoris"
attr(adsl$bl_valvular, "label")     <- "Valvular heart disease"
attr(adsl$bl_card_arrest, "label")  <- "Cardiac arrest"
attr(adsl$bl_AFib, "label")         <- "Atrial fibrillation"
attr(adsl$bl_stroke, "label")       <- "Stroke"
attr(adsl$bl_PAD, "label")          <- "Peripheral artery disease"
attr(adsl$bl_Hypertension, "label") <- "Hypertension"
attr(adsl$bl_Diabetes, "label")     <- "Diabetes mellitus"
attr(adsl$bl_Conduct_dis, "label")  <- "Cardiac conduction disorder"
attr(adsl$bl_VTE, "label")          <- "Venous thromboembolism"






KPI_2009_2022 <- data.frame(
  lisayr = 2009:2022,
  KPI = c(
    299.66, 303.46, 311.43, 314.20, 314.06,
    313.49, 313.35, 316.43, 322.11, 328.40,
    334.26, 335.92, 343.19, 371.91
  )
) %>%
  mutate(KPI_factor = 299.66 / KPI)




mutate(
  across(all_of(paste0("OPD", 1:30)), ~ as.Date(.x, format = "%Y-%m-%d")),
  EVENTDATE = coalesce(
    do.call(pmin, c(pick(all_of(paste0("OPD", 1:30))), na.rm = TRUE)),
    UTDATUM
  )
)





data_check <- data %>%
  group_by(LopNr) %>%
  mutate(
    match = any(
      sapply(seq_along(UTDATUM), function(i) {
        any(UTDATUM[i] == INDATUM[-i], na.rm = TRUE)
      })
    )
  ) %>%
  ungroup()


CreateEVENTS(
  ICD10     = "^(309[4678]|308[123]|3161|3157|317[01]|FPE(00|10|20|26|40|96)|FPF(00|10|20|96)|FPG(10|30|33|36|50|96)|DF013|ZXG[45]0)",
  DISEASE   = "Pacemaker",
  POSITION  = 2,
  YEAR_MIN  = 1996,
  SAVE_PATH = SAVE_PATH
)

ICD10 = "^(FPE(00|10|20|26|40|96)|FPF(00|10|20|96)|FPG(10|30|33|36|50|96)|DF013|ZXG[45]0)"
ICD10 = "^(FQA(00|10|20|30|40|96)|FX[456]0|FXN00)"
ICD10 = "^(DF0(12|17|25|28)|ZXG05)"

scale_color_manual(values = c(
  "Controls"                 = "#F8766D",
  "Renal tumors"             = "#A3A500",
  "Malignant bone tumors"    = "#00BF7D",
  "Radiation to the heart"   = "#00B0F6",
  "Radiation to the abdomen" = "#E69F00"
)) +


compute_ARR_adj <- function(data, outcome, exposure, adjustments,
                            bootstrap = FALSE) {

  data <- data %>%
    mutate(
      tde = .data[[paste0(exposure, "_td")]],
      event = .data[[paste0(outcome, "_event")]]
    )

  # Overall rate ratio: cases vs controls
  fit_overall <- glm(
    event ~ group,
    offset = log(person_time),
    family = poisson(link = "log"),
    data = data
  )

  RR_survivors <- exp(coef(fit_overall)["groupcases"])

  # Function to calculate ARR for bootstrap
  calculate_ARR <- function(dat) {

    # Overall rate ratio
    fit_overall_b <- glm(
      event ~ group,
      offset = log(person_time),
      family = poisson(link = "log"),
      data = dat
    )

    RR_survivors_b <- exp(coef(fit_overall_b)["groupcases"])

    # Adjusted rate ratio for cases
    dat_cases <- filter(dat, group == "cases")

    covs <- adjustments[[outcome]][["cases"]]

    formula <- as.formula(
      paste("event ~", paste(c("tde", covs), collapse = " + "))
    )

    fit_b <- glm(
      formula,
      offset = log(person_time),
      family = poisson(link = "log"),
      data = dat_cases
    )

    RR_b <- exp(coef(fit_b)["tde"])

    # Proportion of events among exposed
    pd_b <- mean(
      dat_cases$tde[dat_cases$event == 1],
      na.rm = TRUE
    )

    # Population attributable fraction
    PAF_b <- pd_b * (RR_b - 1) / RR_b

    # Attributable relative rate
    as.numeric(PAF_b * RR_survivors_b * 100)
  }

  # Bootstrap ARR confidence interval
  ARR_CI <- c(NA_real_, NA_real_)

  if (bootstrap) {

    # File path for exposure-specific CI
    ci_path <- file.path(
      "/safe/data/APNC/SCARLET Riskfactors/Data/Results",
      paste0(exposure, "_CI.xlsx")
    )

    if (file.exists(ci_path)) {

      # Load previously calculated CI
      ARR_CI <- as.numeric(
        readxl::read_xlsx(ci_path)$ARR_CI
      )

      if (length(ARR_CI) != 2 || any(!is.finite(ARR_CI))) {
        stop("Saved ARR confidence interval is invalid: ", ci_path)
      }

    } else {

      # Split data once by individual
      if (!"LopNr" %in% names(data) || anyNA(data$LopNr)) {
        stop("LopNr must exist and contain no missing values.")
      }

      split_data <- split(data, data$LopNr)
      ids <- names(split_data)

      # Bootstrap with 5000 samples
      ARR_boot <- replicate(5000, {

        sampled_ids <- sample(
          ids,
          length(ids),
          replace = TRUE
        )

        boot_data <- do.call(
          rbind,
          unname(split_data[sampled_ids])
        )

        tryCatch(
          calculate_ARR(boot_data),
          error = function(e) NA_real_
        )
      })

      # Keep valid bootstrap estimates
      valid_boot <- ARR_boot[is.finite(ARR_boot)]

      if (length(valid_boot) < 5000) {
        warning(
          length(valid_boot),
          " of 5000 bootstrap estimates were valid."
        )
      }

      if (length(valid_boot) < 5000) {
        stop(
          "Not all bootstrap iterations produced valid estimates. ",
          "Confidence interval was not saved."
        )
      }

      # Calculate 95% percentile confidence interval
      ARR_CI <- as.numeric(
        quantile(
          valid_boot,
          probs = c(0.025, 0.975),
          na.rm = TRUE
        )
      )

      # Save confidence interval
      dir.create(dirname(ci_path), recursive = TRUE, showWarnings = FALSE)

      writexl::write_xlsx(
        x = tibble(
          ARR_CI = ARR_CI
        ),
        path = ci_path
      )
    }
  }

  # Calculate adjusted estimates for cases and controls
  map_dfr(c("cases", "controls"), function(grp) {

    dat <- filter(data, group == grp)

    # Select adjustment variables for outcome and group
    covs <- adjustments[[outcome]][[grp]]

    # Create model formula
    formula <- as.formula(
      paste("event ~", paste(c("tde", covs), collapse = " + "))
    )

    # Poisson regression
    fit <- glm(
      formula,
      offset = log(person_time),
      family = poisson(link = "log"),
      data = dat
    )

    # Rate ratio and 95% CI
    RR <- exp(coef(fit)["tde"])
    CI <- exp(confint.default(fit)["tde", ])

    # Proportion of events among exposed
    pd <- mean(
      dat$tde[dat$event == 1],
      na.rm = TRUE
    )

    # Population attributable fraction
    PAF <- pd * (RR - 1) / RR

    # Attributable relative rate (cases only)
    ARR <- if (grp == "cases") {
      PAF * RR_survivors * 100
    } else {
      NA_real_
    }

    # Format ARR with optional bootstrap CI
    ARR_result <- if (grp == "cases") {

      if (bootstrap) {
        sprintf(
          "%.2f (%.2f-%.2f)",
          ARR, ARR_CI[1], ARR_CI[2]
        )
      } else {
        sprintf("%.2f", ARR)
      }

    } else {
      ""
    }

    tibble(
      Outcome = outcome,
      Exposure = exposure,
      Group = ifelse(grp == "cases", "Cases", "Controls"),
      RR = sprintf("%.2f (%.2f-%.2f)", RR, CI[1], CI[2]),
      pd = sprintf("%.2f", pd * 100),
      PAF = sprintf("%.2f", as.numeric(PAF) * 100),
      ARR = ARR_result
    )
  })
}



















compute_ARR_adj <- function(data, outcome, exposure, adjustments,
                            bootstrap = FALSE) {

  data <- data %>%
    mutate(
      tde = .data[[paste0(exposure, "_td")]],
      event = .data[[paste0(outcome, "_event")]]
    )

  # Function to calculate ARR
  calculate_ARR <- function(dat) {

    # Overall rate ratio: cases vs controls
    fit_overall <- glm(
      event ~ group,
      offset = log(person_time),
      family = poisson(link = "log"),
      data = dat
    )

    RR_survivors <- exp(coef(fit_overall)["groupcases"])

    # Adjusted rate ratio for cases
    dat_cases <- filter(dat, group == "cases")
    covs <- adjustments[[outcome]][["cases"]]

    formula <- reformulate(
      c("tde", covs),
      response = "event"
    )

    fit <- glm(
      formula,
      offset = log(person_time),
      family = poisson(link = "log"),
      data = dat_cases
    )

    RR <- exp(coef(fit)["tde"])

    pd <- mean(dat_cases$tde[dat_cases$event == 1], na.rm = TRUE)

    PAF <- pd * (RR - 1) / RR

    as.numeric(PAF * RR_survivors * 100)
  }

  # Bootstrap ARR confidence interval
  ARR_CI <- c(NA_real_, NA_real_)

  if (bootstrap) {

    split_data <- split(data, data$LopNr)
    ids <- names(split_data)

    ARR_boot <- replicate(5000, {

      sampled_ids <- sample(ids, length(ids), replace = TRUE)

      boot_data <- do.call(
        rbind,
        unname(split_data[sampled_ids])
      )

      tryCatch(
        calculate_ARR(boot_data),
        error = function(e) NA_real_
      )
    })

    valid_boot <- ARR_boot[is.finite(ARR_boot)]

    if (length(valid_boot) < 5000) {
      warning(
        length(valid_boot),
        " of 5000 bootstrap estimates were valid."
      )
    }

    if (length(valid_boot) > 0) {
      ARR_CI <- quantile(
        valid_boot,
        probs = c(0.025, 0.975),
        na.rm = TRUE
      )
    }
  }

  # Overall rate ratio: cases vs controls
  fit_overall <- glm(
    event ~ group,
    offset = log(person_time),
    family = poisson(link = "log"),
    data = data
  )

  RR_survivors <- exp(coef(fit_overall)["groupcases"])

  map_dfr(c("cases", "controls"), function(grp) {

    dat <- filter(data, group == grp)

    # Select adjustment variables
    covs <- adjustments[[outcome]][[grp]]

    # Create model formula
    formula <- reformulate(
      c("tde", covs),
      response = "event"
    )

    # Poisson regression
    fit <- glm(
      formula,
      offset = log(person_time),
      family = poisson(link = "log"),
      data = dat
    )

    # Rate ratio and 95% CI
    RR <- exp(coef(fit)["tde"])
    CI <- exp(confint.default(fit)["tde", ])

    # Proportion of events among exposed
    pd <- mean(dat$tde[dat$event == 1], na.rm = TRUE)

    # Population attributable fraction
    PAF <- pd * (RR - 1) / RR

    # Attributable relative rate (cases only)
    ARR <- if (grp == "cases") {
      PAF * RR_survivors * 100
    } else {
      NA_real_
    }

    # Format ARR
    ARR_result <- if (grp == "cases") {

      if (bootstrap) {
        sprintf(
          "%.2f (%.2f-%.2f)",
          ARR, ARR_CI[1], ARR_CI[2]
        )
      } else {
        sprintf("%.2f", ARR)
      }

    } else {
      ""
    }

    tibble(
      Outcome = outcome,
      Exposure = exposure,
      Group = ifelse(grp == "cases", "Cases", "Controls"),
      RR = sprintf("%.2f (%.2f-%.2f)", RR, CI[1], CI[2]),
      pd = sprintf("%.2f", pd * 100),
      PAF = sprintf("%.2f", as.numeric(PAF) * 100),
      ARR = ARR_result
    )
  })
}


















# Define adjustments for each outcome and group
adjustments <- list(
  MACE = list(
    cases = c("age", "sex"),
    controls = c("age", "sex")
  ),
  Death = list(
    cases = c("age", "sex"),
    controls = c("age", "sex")
  )
)

# Function
calculate_ARR2 <- function(data, outcome, exposure, adjustments) {

  data <- data %>%
    mutate(
      tde = .data[[paste0(exposure, "_td")]],
      event = .data[[paste0(outcome, "_event")]]
    )

  map_dfr(c("cases", "controls"), function(grp) {

    dat <- filter(data, group == grp)

    # Select adjustment variables for outcome and group
    covs <- adjustments[[outcome]][[grp]]

    # Create model formula
    formula <- reformulate(
      c("tde", covs),
      response = "event"
    )

    # Poisson regression
    fit <- glm(
      formula,
      offset = log(person_time),
      family = poisson(link = "log"),
      data = dat
    )

    # Rate ratio and 95% CI
    RR <- exp(coef(fit)["tde"])
    CI <- exp(confint.default(fit)["tde", ])

    # Proportion of events among exposed
    pd <- mean(dat$tde[dat$event == 1], na.rm = TRUE)

    # Population attributable fraction
    PAF <- pd * (RR - 1) / RR

    tibble(
      Outcome = outcome,
      Exposure = exposure,
      Group = ifelse(grp == "cases", "Cases", "Controls"),
      RR = sprintf("%.2f (%.2f-%.2f)", RR, CI[1], CI[2]),
      pd = sprintf("%.2f", pd * 100),
      PAF = sprintf("%.2f", as.numeric(PAF) * 100)
    )
  })
}

# Run analyses
result <- calculate_ARR2(
  data = data,
  outcome = "MACE",
  exposure = "Hypertension",
  adjustments = adjustments
)










adsl %>%
  summarise(
    Heart = sum(radiation_heart == 1, na.rm = TRUE),
    TBI = sum(radiation_totalbody == 1, na.rm = TRUE),
    Both = sum(radiation_heart == 1 & radiation_totalbody == 1, na.rm = TRUE),
    Either = sum(radiation_heart == 1 | radiation_totalbody == 1, na.rm = TRUE)
  )





####################################################################################################
#                                                                                                  #
# ¤Name:         Figures2.R                                                                        #
#                                                                                                  #
# ¤Project:      HelenSjöland/Telemedicine                                                         #
#                                                                                                  #
# ¤Description:  Re-create Figure 4 from Hamah Saed                                                #
#                                                                                                  #
# ¤Date:         2026-02-10                                                                        #
#                                                                                                  #
# ¤Author:       Mathis Wiesner                                                                    #
#                                                                                                  #
# ¤Indata:                                                                                         #
#                                                                                                  #
# ¤Outdata:                                                                                        #
#                                                                                                  #
# ¤Output:                                                                                         #
#                                                                                                  #
####################################################################################################
#                                                                                                  #
# ¤Modifications                                                                                   #
#                                                                                                  #
# ¤ModDate     ¤ModInit ¤ModDescription                                                            #
#                                                                                                  #
####################################################################################################
# ¤Library installation
install.packages("readxl")
install.packages("TTR")
install.packages("writexl")
####################################################################################################
# ¤Library 
library(readxl)
library(ggplot2)
library(dplyr)
library(writexl)
library(haven)
library(tidyr)
library(grid)
library(ggbrace)
library(emmeans)
library(cowplot)
library(broom)
library(sandwich)
library(scales)
####################################################################################################
# ¤Formats
pval_func <- function(pval) {
  case_when(
    is.na(pval)    ~ "",
    pval > 0.1     ~ sprintf("%.2f", pval),
    pval >= 0.0001 ~ sprintf("%.3f", pval),
    TRUE           ~ "<.0001"
  )
}
####################################################################################################
# ¤Read in the data
projroot <- "S:/Cons/HelenSjöland/Telemedicine"
dataroot <- "Q:/Cons/HelenSjöland/Telemedicine"

fullpath1 <- file.path(dataroot, "data", "derived", "adsl_iptw.sas7bdat")
adsl <- read_sas(fullpath1) %>% filter(excl==0) 
adsl <- adsl %>% rename(ATTWGT = `_ATTWGT_`) %>%
                 mutate(trtgrp = relevel(factor(trtgrp), ref = "Controls")) 



#######################################################################
### Implementation of GDMT (any dose) across all four drug classes. ###
#######################################################################

# Define timepoints
timepoints <- c(0, 90, 180, 365, 730)
# Create longitudinal dataset
gdmt_long <- bind_rows(lapply(timepoints, function(day) {
  
  vars <- paste0(c("sglt2h_Target_","MRA_Target_","betablocker_Target_","RAASblock_Target_"),day, "days")
  
  adsl %>%mutate(time = day,
      
      # Alive at this timepoint
      alive = if_all(
        all_of(vars),
        ~ !is.na(.) & . != "Death"
      ),
      
      # All four GDMT classes achieved
      GDMT_target = as.integer(
        if_any(
          all_of(vars),
          ~ . == "gdtm"
        )
      )
    ) %>%
    filter(alive)
}))

#######################################################################
###   IPTW-weighted proportions at index, 12 weeks, 6 - 24 months   ###
#######################################################################
gdmt_percentages <- gdmt_long %>%
  group_by(trtgrp, time) %>%
  summarise(
    N = n(),
    N_GDMT = sum(GDMT_target),
    Percentage = 100 * weighted.mean(
      GDMT_target,
      ATTWGT,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

#######################################################################
###         ORs with 95% CIs are derived from IPTW-weighted         ###
#######################################################################
gdmt_OR <- gdmt_long %>%
  group_by(time) %>%
  group_modify(~ {
    
    fit <- glm(
      GDMT_target ~ trtgrp,
      family = quasibinomial(link = "logit"),
      weights = ATTWGT,
      data = .x
    )
    
    # Robust sandwich variance
    V <- vcovHC(fit, type = "HC0")
    
    tidy(fit) %>%
      filter(term != "(Intercept)") %>%
      mutate(
        SE = sqrt(diag(V)[term]),
        OR = exp(estimate),
        Lower = exp(estimate - 1.96 * SE),
        Upper = exp(estimate + 1.96 * SE),
        p_value = 2 * pnorm(-abs(estimate / SE))
      ) %>%
      select(OR, Lower, Upper, p_value)
    
  }) %>%
  ungroup()


#######################################################################
###                        Create Figure                            ###
#######################################################################

# 1. Prepare percentages
plot_data <- gdmt_percentages %>%
  mutate(
    time = factor(time,levels = c(0, 90, 180, 365, 730),labels = c("Index", "Day 90", "Day 180", "Day 365", "Day 730")),
    trtgrp = factor(trtgrp,levels = c("Controls", "Cases"))
  )

# 2. Calculate p-values from OR and 95% CI
or_data <- gdmt_OR %>%
  mutate(
    time = factor(
      time,
      levels = c(0, 90, 180, 365, 730),
      labels = c("Index", "Day 90", "Day 180", "Day 365", "Day 730")
    ),
    label = sprintf(
      "OR: %.2f (%.2f-%.2f)\nP: %s",
      OR, Lower, Upper, pval_func(p_value)
    )
  )

# 3. Create figure
p <- ggplot(
  plot_data,
  aes(x = time, y = Percentage, fill = trtgrp)
) +
  geom_col(
    position = position_dodge(width = 0.85),
    width = 0.80
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", Percentage)),
    position = position_dodge(width = 0.85),
    vjust = -0.5,
    size = 4,
    fontface = "bold"
  ) +
  geom_text(
    data = or_data,
    aes(x = time, y = 113, label = label),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "italic",
    lineheight = 1.3
  ) +
  scale_fill_manual(
    values = c(
      "Controls" = "#D95F02",
      "Cases" = "#0072B2"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 120),
    breaks = seq(0, 120, 30),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Timepoint",
    y = "Percentage",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey90"),
    axis.text.x = element_text(
      face = "bold",
      color = "grey30"
    ),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(15, 15, 10, 15)
  )

p

ggsave(
  "S:/Cons/HelenSjöland/Telemedicine/Output/figure2.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 600,
  device = ragg::agg_png
)

