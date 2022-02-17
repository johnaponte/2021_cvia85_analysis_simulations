#' ---
#' title: CVIA-085 Analysis of precision simulations using t-test
#' author: John Aponte
#' date: "2021-01-05"
#' ---

#'
#' ## Introduction
#'
#' This program make the analysis of the simulations with different precision
#' using the weighted t-test using as weights the inverse of the variance proposed
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
#' ## Setup
#'
#+ echo = T, result = "hide", message = F, warning = F
starttime_ <- Sys.time()
library(DBI)
library(repana)
library(plyr)
library(tidyverse)
library(broom)
library(weights)
options(knitr.duplicate.label = "allow")
con <- get_con()


#'
#' ## Analysis
#'
#+ echo = T

# Load the simulation matrix
sim_matrix_prec <- dbReadTable(con,"sim_matrix_prec")


# Unweighted t-test ####
utt_prec <-
  sim_matrix_prec %>%
  ddply(
    .(idsim),
    function(x){
      cat("u idsim: ",x$idsim,"\n")
      tblname <- paste0("SP_", formatC(x$idsim,width = 4, format = "d", flag = "0"))
      dba <- dbReadTable(con, tblname)
      ana_summary <-
        ddply(
          dba,
          .(idtrial),
          function(y){
            cat("u", x$idsim, y$idtrial[1],"\n")
            plist <-
              list(
                #d_a0_d1-d_a0_d0
                "d_a0_d" = (y$inf_a0_d1/y$mosq_a0_d1) - (y$inf_a0_d0/y$mosq_a0_d0),
                #d_a1_d1-d_a1_d0
                "d_a1_d" = (y$inf_a1_d1/y$mosq_a1_d1) - (y$inf_a1_d0/y$mosq_a1_d0),
                #d_a1_d0-d_a0_d0
                "d_d0_a" = (y$inf_a1_d0/y$mosq_a1_d0) - (y$inf_a0_d0/y$mosq_a0_d0),
                #d_a1_d1-d_a0_d1
                "d_d1_a" = (y$inf_a1_d1/y$mosq_a1_d1) - (y$inf_a0_d1/y$mosq_a0_d1)
              )
            res <-
              ldply(
                plist,
                function(z){
                  z <- z[is.finite(z)]

                  t.test(z, mu = 0) %>% tidy()
                }
              )
            res
          })
      ana_summary
    })

update_table(con, utt_prec, "Unweighted ttest")

# Weighted t-test ####

tidy_wgt.t <- function(xx){
  # tidy function to convert the results in a tibble
  tibble(
    estimate = xx$additional["Difference"],
    statistics = xx$coefficients["t.value"],
    p.value = xx$coefficients["p.value"],
    parameter = xx$coefficients["df"] ,
    conf.low = xx$additional["Difference"] - qt(0.975,df = xx$coefficients["df"])*xx$additional["Std. Err"],
    conf_high = xx$additional["Difference"] + qt(0.975,df = xx$coefficients["df"])*xx$additional["Std. Err"],
    method = xx$test,
    alternative = "two.sided"
  )}

## Analysis of weighted t
wtt_prec <-
  sim_matrix_prec %>%
  ddply(
    .(idsim),
    function(x){
      cat("w idsim: ",x$idsim,"\n")
      tblname <- paste0("SP_", formatC(x$idsim,width = 4, format = "d", flag = "0"))
      dba <- dbReadTable(con, tblname)
      ana_summary <-
        ddply(
          dba,
          .(idtrial),
          function(y){
            cat("w", x$idsim, y$idtrial[1],"\n")
            pa0d0 <- (y$inf_a0_d0+1)/(y$mosq_a0_d0+2)
            pa0d1 <- (y$inf_a0_d1+1)/(y$mosq_a0_d1+2)
            pa1d0 <- (y$inf_a1_d0+1)/(y$mosq_a1_d0+2)
            pa1d1 <- (y$inf_a1_d1+1)/(y$mosq_a1_d1+2)
            d_a0_d <- pa0d1-pa0d0
            w_a0_d <- 1/((pa0d1*(1-pa0d1))/(y$mosq_a0_d1+2) + (pa0d0*(1-pa0d0))/(y$mosq_a0_d0+2))

            d_a1_d <- pa1d1-pa1d0
            w_a1_d <- 1/((pa1d1*(1-pa1d1))/(y$mosq_a1_d1+2) + (pa1d0*(1-pa1d0))/(y$mosq_a1_d0+2))

            d_d0_a <- pa1d0-pa0d0
            w_d0_a <- 1/((pa1d0*(1-pa1d0))/(y$mosq_a1_d0+2) + (pa0d0*(1-pa0d0))/(y$mosq_a0_d0+2))

            d_d1_a <- pa1d1-pa0d1
            w_d1_a <- 1/((pa1d1*(1-pa1d1))/(y$mosq_a1_d1+2) + (pa0d1*(1-pa0d1))/(y$mosq_a0_d1+2))

            res <-
              wtd.t.test(d_a0_d, y=0,weight = w_a0_d) %>% tidy_wgt.t() %>%
              bind_rows(
                wtd.t.test(d_a1_d, y=0,weight = w_a1_d) %>% tidy_wgt.t()
              ) %>%
              bind_rows(
                wtd.t.test(d_d0_a, y=0,weight = w_d0_a) %>% tidy_wgt.t()
              ) %>%
              bind_rows(
                wtd.t.test(d_d1_a, y=0,weight = w_d1_a) %>% tidy_wgt.t()
              ) %>%
              mutate(".id" = c("d_a0_d","d_a1_d","d_d0_a","d_d1_a"))
            res
          })
      ana_summary
    })

update_table(con, wtt_prec, "Weighted ttest")

dbDisconnect(con)
#'
#' ## Session Info
#'
#+  echo = F
cat("Execution date: " , format(starttime_ ),"\n")
cat("Execution time: ", round(difftime(Sys.time(),starttime_,units = "min"),2), "min\n")
print(sessionInfo(), locale = F)
