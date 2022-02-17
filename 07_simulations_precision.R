#' ---
#' title: CVIA-085 simulations for precision
#' author: John Aponte
#' date: "2021-11-10"
#' ---

#'
#' ## Introduction
#'
#' This program simulates different levels of precision, assuming there is
#' no difference between days or between assays
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
con <- get_con()
# if true run the simulation
simulate = TRUE

# helper functions and values
to_mu <- function(shape1,shape2){shape1/(shape1+shape2)}
to_phi <- function(shape1,shape2){shape1 + shape2}
to_shape1 <- function(mu,phi){mu*phi}
to_shape2 <- function(mu,phi){phi*(1-mu)}
to_odds <- function(p){p/(1-p)}
to_p <- function(odds){odds/(1+odds)}

to_icc <- function(shape1,shape2){1/(shape1+shape2+1)}
to_shape1_icc <- function(mu,icc){mu*(1-icc)/icc}
to_shape2_icc <- function(mu,icc){(1-mu)*(1-icc)/icc}
to_icc_phi <- function(phi){1/(phi + 1)}
to_phi_icc <- function(icc){1/icc - 1}

#'
#' ## Analysis
#'
#+ echo = T
#Simulate one trial
#
# @param mu mean probability for mosquito infection (assay0 day0)
# @param phi precision parameter for the mosquito infection
# @param nsubjects number of subjects to be simulated
# @param nmosquitoes number of mosquitoes that start assays
# @param mosq_mort_low lower limit of mosquito mortality
# @param mosq_mort_upper upper limit of mosquito mortality
# @param day_OR Odds Ratio of $day_1/day_2$
# @param assay_OR Odds Ratio of $assay_1/assay_0$
# @param precision ratio of $assay_1/assay_0$
# @param day_corr Correlation between $day_1$ and $day_0$
# @import copula
# @import magrittr
# @import tibble
# @return a tibble with the number of positive and total surviving mosquitoes by
# assay and day, with one row per subject
simul_one_prec <- function(
  mu,
  phi,
  nsubjects,
  nmosquitoes,
  mosq_mort_low,
  mosq_mort_high,
  day_OR,
  assay_OR,
  precision_ratio,
  day_corr
){

  # Step 1: Simulate the prior mosquito survival probability
  priorp <- runif(nsubjects, min = 1-mosq_mort_high, max = 1-mosq_mort_low)

  # Step 2: Simulate the surviving mosquitoes in a matrix nsubjects x 4
  mosq <- matrix(
    rbinom(nsubjects*4, nmosquitoes, priorp),
     ncol = 4,
     byrow = F,
     dimnames = list(1:nsubjects,c("mosq_a0_d0", "mosq_a1_d0", "mosq_a0_d1", "mosq_a1_d1")))

  # Step 3: Estimate mu parameters for the different assays
  a0_d0 <- mu
  a0_d1 <- to_p(to_odds(mu)*day_OR)
  a1_d0 <- to_p(to_odds(mu)*assay_OR)
  a1_d1 <- to_p(to_odds(mu)*day_OR *assay_OR)

  phi_a0 <- phi
  phi_a1 <- phi*precision_ratio

  # Step 4: Create the copula between days according to the correlation
  # In R, beta distribution is simulated in shape1, shape2 parametrisation
  betacopula <- mvdc(
    copula = normalCopula(
      param = c(1, day_corr, day_corr, day_corr,day_corr,1),
      dim = 4,
      dispstr = "un"
    ),
    margins = c("beta", "beta","beta","beta"),
    paramMargins = list(
      list(shape1 = to_shape1(a0_d0,phi_a0), shape2 = to_shape2(a0_d0,phi_a0)),
      list(shape1 = to_shape1(a1_d0,phi_a1), shape2 = to_shape2(a1_d0,phi_a1)),
      list(shape1 = to_shape1(a0_d1,phi_a0), shape2 = to_shape2(a0_d1,phi_a0)),
      list(shape1 = to_shape1(a1_d1,phi_a1), shape2 = to_shape2(a1_d1,phi_a1))
    )
  )

  #' Step 5: Simulate the probability of infection by subject and assay in a matrix
  pinfection <- rMvdc(nsubjects, betacopula)
  infected <- matrix(
    rbinom(nsubjects*4, mosq,pinfection),
    ncol = 4,
    byrow = F,
    dimnames = list(1:nsubjects,c("inf_a0_d0", "inf_a1_d0", "inf_a0_d1", "inf_a1_d1")))

  #' Step 6: Return a tibble with the simulated infected and surviving mosquitoes
  as_tibble(cbind(infected,mosq))
 }

# Define the simulation matrix
sim_matrix_prec <-
  expand.grid(
    mu = 0.1685,
    phi = 5.1535,
    nsubjects = c(30,1000),
    nmosquitoes =30,
    mosq_mort_high = .80,
    mosq_mort_low = .20,
    day_OR = seq(1),
    assay_OR = seq(1),
    precision_ratio = c(1,1.3,1.5,2),
    day_corr = c(0,0.5,1),
    nsimul = 10000
  ) %>%
  mutate(idsim = 1:n())

update_table(con,sim_matrix_prec,format(date()))

if (simulate) {
# Perform the simulation
set.seed(23456)
sim_matrix_prec %>%
  d_ply(
    .(idsim),
    .progress = "text",
    function(x){
      sdf <- adply(
        c(1:x$nsimul),
        .id = "idtrial",
        .margins = 1,
        function(y){
          simul_one_prec(
            mu = x$mu,
            phi = x$phi,
            nsubjects = x$nsubject,
            nmosquitoes = x$nmosquitoes,
            mosq_mort_high = x$mosq_mort_high,
            mosq_mort_low = x$mosq_mort_low,
            day_OR = x$day_OR,
            assay_OR = x$assay_OR,
            precision_ratio = x$precision_ratio,
            day_corr = x$day_corr
          )
        }
      )
      tblname <- paste0("SP_", formatC(x$idsim,width = 4, format = "d", flag = "0"))
      update_table(con, sdf, format(date()),tablename = tblname)
    }
  )

}
dbDisconnect(con)

#'
#' ## Session Info
#'
#+  echo = F
cat("Execution date: " , format(starttime_ ),"\n")
cat("Execution time: ", round(difftime(Sys.time(),starttime_,units = "min"),2), "min\n")
print(sessionInfo(), locale = F)
