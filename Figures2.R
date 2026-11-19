

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
