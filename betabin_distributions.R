#' ---
#' title: Betabinomial distributions with different dispersion values
#' author: John Aponte
#' date: "2022-01-07"
#' ---

#'
#' ## Introduction
#'
#' This program make the density graphs for beta-binomial distributions
#' with different dispersion values for the same mean proportion.
#' The previous epidemiological study estimations for the mean proportion is
#' 16.85% and dispersion 5.1535 (ICC 16.25%)
#'
#' ## Parametrization of the beta-binomial distribution
#'
#' In addition of the number of trials parameter, there are two parameters for the
#' beta part of the distribution that modify the binomial part of the distribution.
#'
#' 1. canonical parametrization based on two shape parameters: $shape_1$ and $shape_2$
#'
#' 2. Based on precision parameter: $\mu$ and $\phi$
#'
#'      $\mu = shape1/(shape_1 + shape_2)$
#'
#'      $\phi = shape_1 + shape_2$
#'
#'      $shape_1 =  \mu\phi$
#'
#'      $shape_2 = (1-\mu)\phi$
#'
#' 3. Based on intra-class correlation parameter: $\mu$ and $\rho$
#'
#'      $\mu = shape1/(shape_1 + shape_2)$
#'
#'      $\rho = 1/(shape_1 + shape_2 + 1)$
#'
#'      $shape_1 = \mu\phi$
#'
#'      $shape_2 = (1-mu)(1-\rho)/\rho$
#'
#'
#' The relation between the precision parameter and the icc is:
#'
#' $\rho = 1/(\phi-1)$
#'
#' The lower the intra-class correlation or higher the precision the distribution
#' is more similar to the binomial distribution.
#'
#'
#' ## Setup
#'
#+ echo = T, result = "hide", message = F, warning = F
starttime_ <- Sys.time()
library(DBI)
library(repana)
library(plyr)
library(tidyverse)
library(convdistr)
options(knitr.duplicate.label = "allow")

#'
#' ## Analysis
#'
#+ echo = T
od = 5.1535
dfp <- expand.grid(
  mu = 0.1685,
  od = c(od*0.1,od*0.3, od*0.5,od*0.75, od, od*1.25, od*2.5, od*3, od*10),
  size = 30) %>%
  mutate(icc = 1/(od +1))


simdistr <-
  dfp %>%
  ddply(
    .(mu, od, icc, size),
    function(x){
      rfunc(new_BETABINOMIAL_od(x$size, x$mu, x$od), 1000)
  })

simdistr %>%
  mutate(iccfmt = round(icc*100,1)) %>%
  ggplot() +
  aes(x = rvar) +
  geom_density()+
  facet_wrap(.~iccfmt) +
  ggtitle(
    "Beta-binomial distributions",
    subtitle ="Different ICC values with a mean proportion of 16.8%" )

#'
#' ## Session Info
#'
#+  echo = F
cat("Execution date: " , format(starttime_ ),"\n")
cat("Execution time: ", round(difftime(Sys.time(),starttime_,units = "min"),2), "min\n")
print(sessionInfo(), locale = F)
