
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
