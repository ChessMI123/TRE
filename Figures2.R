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
