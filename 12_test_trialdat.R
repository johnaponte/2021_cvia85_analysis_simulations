#' ---
#' title: Test results with trial data
#' author: John Aponte
#' date: "2022-06-01"
#' ---

#'
#' ## Introduction
#'
#' To compare SAS and R beta binomial models
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
library(rjags)
library(coda)
library(ggmcmc)
options(knitr.duplicate.label = "allow")

axl05 <- read.csv("_data/axl05.csv")


model <- "
data {
  Nsubj = length(subject)
}

model {
  # ################################################################
  # BetaBinomial model for rjags
  # logistic model for the mean probability
  # loglinear model for the precision parameter
  ##################################################################

  # Parameters that modify the probability in a logit model
  alpha ~ dnorm(0, 0.01)
  mu00 <- ilogit(alpha)

  # Parameters that modify the dispersion in a log linear model
  gamma ~ dnorm(0, 0.01)

  phi00 <- exp(gamma)

  # Conversion from mu and phi to shape_1 and shape_2
  shape1_00 <- mu00*phi00
  shape2_00 <- (1-mu00)*phi00

  for (i in 1:Nsubj){

    # Probabilities for the binomial model
    prob_00[i] ~ dbeta(shape1_00, shape2_00)

    # Binomial model
    pos_00[i] ~ dbinom(prob_00[i], mosq_00[i])

  }
}
"

# Extract the data for baseline dsfa
data_wider <-
  axl05 %>%
  filter(AVISIT == "BASELINE" & LBANMETH == "DMFA" & NMOC > 0) %>%
  mutate(subject = seq(1,n()), pos_00 =NPOC, mosq_00 = NMOC) %>%
  select(subject, pos_00,mosq_00)

# Prepare the model
betabinom_jags <-
  jags.model(
    textConnection(model),
    data = data_wider,
    n.chains = 4,
    n.adapt = 1000,
    inits = list(
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 1111),
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 2222),
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 3333),
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 4444)
    )
  )

# Data example
set.seed(123456)

# Simulate the posterior
latent_sim <-
  coda.samples(
    model = betabinom_jags,
    variable.names = c("shape1_00","shape2_00", "mu00", "phi00" ),
    n.iter =  100000,
    thin = 1)

#+ echo = T
# Extract and Analyze the posterior
latent_sum <-  summary(latent_sim)
eff_size <-  effectiveSize(latent_sim)
latent_chains <- ggs(latent_sim) %>%
  filter(Parameter %in% c("shape1_00","shape2_00", "mu00", "phi00"))

#' ## Results
res <- cbind(latent_sum[[1]], latent_sum[[2]], eff_size)
res[c("shape1_00","shape2_00", "mu00", "phi00"), c("eff_size", "Mean", "SD","50%","2.5%","97.5%")]

#' ## Revision of the convergence of the model
#+ out.width = "100%", dpi = 300 , fig.width = 6, fig.height = 8
ggmcmc(latent_chains, param_page = 6, file = NULL)





#'
#' ## Session Info
#'
#+  echo = F
cat("Execution date: " , format(starttime_ ),"\n")
cat("Execution time: ", round(difftime(Sys.time(),starttime_,units = "min"),2), "min\n")
print(sessionInfo(), locale = F)
