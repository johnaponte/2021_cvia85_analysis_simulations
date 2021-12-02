#' ---
#' title: Summary of the analysis of the simulations
#' author: John Aponte
#' date: "2021-11-16"
#' ---

#'
#' ## Introduction
#' This program summary the results of the t-test and weighted t-test analysis
#' of the simulations made for CVIA-085

#'
#' ## Setup
#'
#+ echo = T, result = "hide", message = F, warning = F
starttime_ <- Sys.time()
library(DBI)
library(repana)
library(plyr)
library(tidyverse)
library(gridExtra)
options(knitr.duplicate.label = "allow")
con <- get_con()
#'
#' ## Analysis
#'
#+ echo = T
sim_matrix <- dbReadTable(con, "sim_matrix")
utt <- dbReadTable(con, "utt")
wtt <- dbReadTable(con, "wtt")

# No weighted t-test ####
power <- utt %>%
  mutate(sig = ifelse(p.value < 0.05, 1, 0)) %>%
  filter(is.finite(sig)) %>%
  rename(id = .id) %>%
  group_by(idsim, id) %>%
  summarise(power = mean(sig, na.rm = T)) %>%
  left_join (sim_matrix, by = "idsim")


# Difference between assay 1 and assay 0 at day 0
d_d0_a <- power %>%
  filter(id == "d_d0_a") %>%
  ggplot() +
  aes(
    x = assay_OR,
    y = power,
    linetype = as.factor(day_corr),
    color = as.factor(day_OR)
  ) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = 2, color = "gray") +
  facet_wrap( ~ nsubjects) +
  scale_x_continuous("Assay OR") +
  scale_y_continuous("Empirical power", breaks = c(0,0.05,0.20,0.40,0.60,0.80,1)) +
  scale_color_discrete("Day OR") +
  scale_linetype_discrete("Correlation\nbetween days") +
  ggtitle("Power to detect a difference between Assays at Day 0")

# Difference between assay1 and assay 0 at day 1
d_d1_a <- power %>%
  filter(id == "d_d1_a") %>%
  ggplot() +
  aes(
    x = assay_OR,
    y = power,
    linetype = as.factor(day_corr),
    color = as.factor(day_OR)
  ) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = 2, color = "gray") +
  facet_wrap( ~ nsubjects) +
  scale_x_continuous("Assay OR") +
  scale_y_continuous("Empirical power", breaks = c(0,0.05,0.20,0.40,0.60,0.80,1)) +
  scale_color_discrete("Day OR") +
  scale_linetype_discrete("Correlation\nbetween days") +
  ggtitle("Power to detect a difference between Assays at Day 1")

# Difference betwen day 1 and day 0 for assay 0
d_a0_d <- power %>%
  filter(id == "d_a0_d") %>%
  ggplot() +
  aes(
    x = day_OR,
    y = power,
    linetype = as.factor(day_corr),
    color = as.factor(assay_OR)
  ) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = 2, color = "gray") +
  facet_wrap( ~ nsubjects) +
  scale_x_continuous("Day OR") +
  scale_y_continuous("Empirical power", breaks = c(0,0.05,0.20,0.40,0.60,0.80,1)) +
  scale_color_discrete("Assay OR") +
  scale_linetype_discrete("Correlation\nbetween days") +
  ggtitle("Power to detect a difference between days for  assay 0")

# Difference between day 1 and day 0 for assay 1
d_a1_d <- power %>%
  filter(id == "d_a1_d") %>%
  ggplot() +
  aes(
    x = day_OR,
    y = power,
    linetype = as.factor(day_corr),
    color = as.factor(assay_OR)
  ) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = 2, color = "gray") +
  facet_wrap( ~ nsubjects) +
  scale_x_continuous("Day OR") +
  scale_y_continuous("Empirical power", breaks = c(0,0.05,0.20,0.40,0.60,0.80,1)) +
  scale_color_discrete("Assay OR") +
  scale_linetype_discrete("Correlation\nbetween days") +
  ggtitle("Power to detect a difference between days for Assay 1")


plots <- marrangeGrob(
  list(d_a0_d, d_a1_d, d_d0_a, d_d1_a),
  nrow = 2,
  ncol = 2,
  top = "Analysis using T-test"
)


# Weighted test ####
wpower <- wtt %>%
  mutate(sig = ifelse(p.value < 0.05, 1, 0)) %>%
  filter(is.finite(sig)) %>%
  rename(id = .id) %>%
  group_by(idsim, id) %>%
  summarise(power = mean(sig)) %>%
  left_join (sim_matrix, by = "idsim")

# Difference between assay 1 and assay 0 at day 0
wd_d0_a <- wpower %>%
  filter(id == "d_d0_a") %>%
  ggplot() +
  aes(
    x = assay_OR,
    y = power,
    linetype = as.factor(day_corr),
    color = as.factor(day_OR)
  ) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = 2, color = "gray") +
  facet_wrap( ~ nsubjects) +
  scale_x_continuous("Assay OR") +
  scale_y_continuous("Empirical power", breaks = c(0,0.05,0.20,0.40,0.60,0.80,1)) +
  scale_color_discrete("Day OR") +
  scale_linetype_discrete("Correlation\nbetween days") +
  ggtitle("Power to detect a difference between Assays at Day 0")

# Difference between assay1 and assay 0 at day 1
wd_d1_a <- wpower %>%
  filter(id == "d_d1_a") %>%
  ggplot() +
  aes(
    x = assay_OR,
    y = power,
    linetype = as.factor(day_corr),
    color = as.factor(day_OR)
  ) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = 2, color = "gray") +
  facet_wrap( ~ nsubjects) +
  scale_x_continuous("Assay OR") +
  scale_y_continuous("Empirical power", breaks = c(0,0.05,0.20,0.40,0.60,0.80,1)) +
  scale_color_discrete("Day OR") +
  scale_linetype_discrete("Correlation\nbetween days") +
  ggtitle("Power to detect a difference between Assays at Day 1")

# Difference betwen day 1 and day 0 for assay 0
wd_a0_d <- wpower %>%
  filter(id == "d_a0_d") %>%
  ggplot() +
  aes(
    x = day_OR,
    y = power,
    linetype = as.factor(day_corr),
    color = as.factor(assay_OR)
  ) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = 2, color = "gray") +
  facet_wrap( ~ nsubjects) +
  scale_x_continuous("Day OR") +
  scale_y_continuous("Empirical power", breaks = c(0,0.05,0.20,0.40,0.60,0.80,1)) +
  scale_color_discrete("Assay OR") +
  scale_linetype_discrete("Correlation\n between days") +
  ggtitle("Power to detect a difference between days for  assay 0")

# Difference between day 1 and day 0 for assay 1
wd_a1_d <- wpower %>%
  filter(id == "d_a1_d") %>%
  ggplot() +
  aes(
    x = day_OR,
    y = power,
    linetype = as.factor(day_corr),
    color = as.factor(assay_OR)
  ) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = 2, color = "gray") +
  facet_wrap( ~ nsubjects) +
  scale_x_continuous("Day OR") +
  scale_y_continuous("Empirical power", breaks = c(0,0.05,0.20,0.40,0.60,0.80,1)) +
  scale_color_discrete("Assay OR") +
  scale_linetype_discrete("Correlation\nbetween days") +
  ggtitle("Power to detect a difference between days for Assay 1")


wplots <- marrangeGrob(
  list(wd_a0_d, wd_a1_d, wd_d0_a, wd_d1_a),
  nrow = 2,
  ncol = 2,
  top = "Analysis using Weighted T-test"
)
pdf(file = "reports/cvia_085_power_analysis.pdf",
    width = 16,
    height = 11)
plots
wplots
dev.off()

dbDisconnect(con)
#'
#' ## Session Info
#'
#+  echo = F
cat("Execution date: " , format(starttime_), "\n")
cat("Execution time: ", round(difftime(Sys.time(), starttime_, units = "min"), 2), "min\n")
print(sessionInfo(), locale = F)
