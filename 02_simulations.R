#' ---
#' title: CVIA-085 simulations taking into account a lower number of mosquitoes
#' author: John Aponte
#' date: "2021-11-10"
#' ---

#'
#' ## Introduction
#'
#' CVIA-085 is a phase0 trial to evaluate the feasibility of two consecutive DSFA
#' and DMFA assays in subjects PCR positive for gametocytemia. In this assays
#' 30 mosquitoes are included in each assay each day with the expectancy that 25
#' will survive at day 9 to evaluate the presence of Oocysts.
#' The trial wants to evaluate the variability of the Oocyst prevalence between
#' days and between assays. In this trial, the mortality of mosquitoes is higher
#' than expected. This file create a simulation of trials to be analyzed using
#' different techniques, taking into account that mosquito mortality is higher than
#' expected.
#'
#' The simulation is based on a beta-binomial distribution for Oocyst prevalence
#' as observed in a previous epidemiological survey, with mean 16.85% and precision 5.1535
#' (See below the estimations from the previous epidemiological study).
#' Day and Assay type will simulate changes on the mean proportion of infected
#' mosquitoes but not on the precision parameter of the distribution. In addition of the
#' effect by assay and day, the simulations will take into account the correlation between
#' the probabilities of the beta-binomial distribution by day, assuming that the
#' prior probability by day will be the same in both assays but it maybe or not
#' correlated between days.
#'
#' Mosquito mortality will be simulated as a binomial distributions having as hyper-
#' parameter a rectangular distribution between 20% and 80%. Each subject will have
#' the same mosquito mortality probability for the four assays, but the mortality
#' will vary between subjects according to the rectangular distribution.
#'
#' The correlation between days and assays will be simulated using a normal
#' copula with beta-distribution margins and the following correlation structure
#'
#' |                   | $assay_0 day_0$  | $assay_1 day_0$ | $assay_0 day_1$ | $assay_1 day_1$  |
#' |:-----------------:|:----------------:|:---------------:|:---------------:|:----------------:|
#' |  $assay_0 day_0$  |  1               |  1              | $\rho$  | $\rho$ |
#' |  $assay_1 day_0$  |  1               |  1              | $\rho$  | $\rho$ |
#' |  $assay_0 day_1$  | $\rho$           |  $\rho$         | 1       | 1      |
#' |  $assay_0 day_1$  | $\rho$           |  $\rho$         | 1       | 1      |
#'
#' ### Estimations from a previous epidemiological study
#'
#' A beta-binomial bayesian model was fitted to the data provided from a previous
#' epidemiological study on the probability of mosquito infection in subjects
#' PCR positive for gametocytes. The summary of the posterior distribution is
#' is presented below:

##                        Mean       SD  Naive SE Time-series SE
## mean (mu)          0.1685 0.007246 3.623e-05      6.316e-05
## precision (phi)    5.1535 0.367023 1.835e-03      6.488e-03
## shape1             0.8676 0.059888 2.994e-04      1.090e-03
## shape2             4.2859 0.320259 1.601e-03      5.415e-03
##

#' The beta-binomial distribution can be presented as a function of two parameters:
#' mean ($\mu$) and precision ($\phi$) or in terms of $shape_1$ and $shape_2$.
#' The association between the two different representations is:
#'
#' $\mu = shape_1/(shape_1 + shape_2)$
#'
#' $\phi = shape_1 + shape_2$
#'
#' $shape_1 = \mu\phi$
#'
#' $shape_1 = \phi(1-\mu)$
#'
#' ## Setup
#'
#+ echo = T, result = "hide", message = F, warning = F
starttime_ <- Sys.time()
library(DBI)
library(repana)
library(plyr)
library(tidyverse)
library(copula)
options(knitr.duplicate.label = "allow")

#'
#' ## Analysis
#'
#+ echo = T





#'
#' ## Session Info
#'
#+  echo = F
cat("Execution date: " , format(starttime_ ),"\n")
cat("Execution time: ", round(difftime(Sys.time(),starttime_,units = "min"),2), "min\n")
print(sessionInfo(), locale = F)
