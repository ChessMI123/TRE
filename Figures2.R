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

