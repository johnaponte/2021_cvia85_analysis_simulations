#' ---
#' title: Documentation for CVIA-085 simulations
#' author: John Aponte
#' date: "2021-11-18"
#' ---

#' ## Introduction
#'
#' CVIA-085 is a phase0 trial to evaluate the feasibility of two consecutive DSFA
#' and DMFA assays in subjects PCR positive for gametocytemia. In this assays
#' 30 mosquitoes are included in each assay each day with the expectancy that 25
#' will survive at day 9 to evaluate the presence of Oocysts.
#' The trial wants to evaluate the variability of the Oocyst prevalence between
#' days and between assays. In this trial, the mortality of mosquitoes is higher
#' than expected. This project create a simulation of trials to be analyzed using
#' different techniques, taking into account that mosquito mortality is higher than
#' expected.
#'
#' ## Simulations
#'
#' The simulations are based on a beta-binomial distribution for Oocyst prevalence
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
#' $\mu_{ij}$ is the mean value for the assay $i$ and the day $j$.
#'
#' The effect of day on the assay is simulated as:
#'
#' OR day for assay 0: $odds(\mu_{01})/odds(\mu_{00})$
#'
#' OR day for assay 1: $odds(\mu_{11})/odds(\mu_{01})$
#'
#' The effect of assay on a day is simulated as:
#'
#' OR assay for day 0: $odds(\mu_{10})/odds(\mu_{00})$
#'
#' OR assay for day 1: $odds(\mu_{11})/odds(\mu_{01}))$
#'
#' ## Estimations from a previous epidemiological study
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
#'
#' ## Analysis of simulations
#'
#' Each simulated trial is analyzed by comparing the change by day and assay
#' on the proportion of infected mosquitoes using a un-weighted t-test and a
#' weighted t-test using as weights the inverse of the variance proposed
#' by Agresti-Caffo (adjusted Wald test), to test the hypothesis that the
#' average difference of proportions is different from 0.
#'
#' The variance of the difference of proportion $\tilde{p}_1 - \tilde{p}_2$ is calculated
#' as:
#'
#' $\frac{\tilde{p}_1(1-\tilde{p}_1)}{n_1+2} + \frac{\tilde{p}_2(1-\tilde{p}_2)}{n_2+2}$
#'
#'where:
#'
#' $\tilde{p}_1 = \frac{X_i+1}{n_1+2}$ and $\tilde{p}_2 = \frac{Y_i+1}{n_2+2}$
#'
#' $X_i$ and $Y_i$ are the number of positive mosquitoes in group 1 and 2 respectively
#'
#' $n_1$ and $n_2$ are the total number of mosquitoes en group 1 and 2 respectively
#'
#' ## Results
#'
#' The following graphs show the empirical power to detect an average difference
#' in proportion different from 0 with a p-value < 0.05 using different OR
#' for assay, day, correlation between days and number of subjects in the trial
#' (N = 30 or N = 1000). 10000 simulations were used for each combination of
#' parameters
#'
#' See file cvia_085_power_analysis.pdf
