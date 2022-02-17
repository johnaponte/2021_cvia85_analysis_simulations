#' ---
#' title: CVIA-085 Width of the confidence interval with different precision
#' author: John Aponte
#' date: "2022-01-05"
#' ---

#'
#' ## Introduction
#'
#' This program analyzes the width of the confidence intervals under the null
#' hypothesis of no difference in the means but difference in precision,
#' based on the weighted t-test comparison. The analysis is restricted
#' to the assays with 30 subjects.
#'
#' ## Setup
#'
#+ echo = F, result = "hide", message = F, warning = F
starttime_ <- Sys.time()
library(DBI)
library(repana)
library(plyr)
library(tidyverse)
library(broom)
library(weights)
library(knitr)
library(kableExtra)
library(DT)
library(writexl)
options(knitr.duplicate.label = "allow")
knitr::opts_chunk$set(echo = FALSE)

con <- get_con()


#'
#' ## Analysis
#'
#+ echo = T

# Load the analyzed data

wtt_prec <- dbReadTable(con,"wtt_prec") %>%
  left_join(dbReadTable(con,"sim_matrix_prec"), by = "idsim")

day_OR1 <-
  wtt_prec %>%
  filter(day_OR == 1) %>%
  filter(nsubjects==30) %>%
  filter(.id %in% c("d_a0_d","d_a1_d")) %>%
  mutate(assay = ifelse(.id == "d_a0_d","Assay 0", "Assay 1")) %>%
  mutate(halfwidth = (conf_high - conf.low)/2 ) %>%
  mutate(sig = ifelse(p.value < 0.05,1,0)) %>%
  group_by(assay, assay_OR) %>%
  arrange(day_corr) %>%
  mutate(idord = seq(1,n()))

summary_OR1 <-
  day_OR1 %>%
  dplyr::group_by(assay, precision_ratio, day_corr) %>%
  mutate(halfwidth = halfwidth*100) %>%
  summarise(
    min = round(min(halfwidth),1),
    p05 = round(quantile(halfwidth, 0.05),1),
    p25 = round(quantile(halfwidth, 0.25),1),
    p50 = round(quantile(halfwidth, 0.50),1),
    mean = round(mean(halfwidth),1),
    p75 = round(quantile(halfwidth, 0.75),1),
    p95 = round(quantile(halfwidth, 0.95),1),
    max = round(max(halfwidth),1),
    power = round(mean(sig),3)) %>%
  arrange(assay, precision_ratio, day_corr)

pl1<- day_OR1 %>%
  ggplot() +
  aes(x = estimate, color = as.factor(day_corr) ) +
  geom_density() +
  facet_grid(assay~precision_ratio) +
  theme(legend.position= "bottom") +
  scale_y_continuous("Density") +
  scale_x_continuous("Difference in proportion Day 1 - Day 0") +
  scale_color_discrete("Correlation between days")+
  ggtitle("Distribution of the estimated difference", subtitle = "by Assay and Assay OR") +
  theme(legend.position = "bottom")
pl1

pl2<- day_OR1 %>%
  ggplot() +
  aes(y =statistics, ymin = conf.low, ymax = conf_high, x = idord, color = as.factor(day_corr)) +
  geom_errorbar(alpha = 0.1) +
  facet_grid(assay~precision_ratio) +
  scale_y_continuous("Proportion difference and 95%CI", limits = c(-0.30,0.30)) +
  scale_x_continuous("Simulation") +
  scale_color_discrete("Correlation between days")+
  ggtitle("Estimated confidence intervals", subtitle = "by Assay and precision_ratio") +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank())
pl2

#' Summary of the width of the 96% confidence intervals for the proportion
#' difference between day 1 and day 0 under the null hypothesis of no difference
#' between days
#+ results = "asis"
#summary_OR1 %>% kable() %>% kable_styling(full_width = FALSE)
datatable(summary_OR1)

pdf(file = "reports/cvia_085_precision_analysis_2.pdf",
    width = 16,
    height = 11)
pl1
pl2
dev.off()

write_xlsx(summary_OR1, path = "reports/cvia_085_precision_analysis_2.xlsx")
#' ## Discussion
#'
#' The precision of the difference between the proportion of infected mosquitoes
#' at day 1 and day 0 under the null hypothesis of no difference in the mean
#' proportion depends on the
#' correlation between the probability of infection at day 1 and day 0, and the
#' precision of the assay. ..TBD


dbDisconnect(con)
rm(con)
#'
#' ## Session Info
#'
#+  echo = F
cat("Execution date: " , format(starttime_ ),"\n")
cat("Execution time: ", round(difftime(Sys.time(),starttime_,units = "min"),2), "min\n")
print(sessionInfo(), locale = F)
