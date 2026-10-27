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

