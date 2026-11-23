
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
