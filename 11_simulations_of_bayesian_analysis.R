#' ---
#' title: Bayesian model for the CVIA-085 trial with both ICC and Mean change
#' author: John Aponte
#' date: "2022-02-16"
#' ---

#'
#' ## Introduction
#'
#' Let's assume the number of infected mosquitoes $n$ from $m$ exposed mosquitos
#' for the subject $i$ the day $j$ and the assay $k$ follows a $Binomial$
#' distributions as:
#'
#' $$ n_{i,j,k} \sim Binomial(\pi_{i,j,k},m_{i,j,k})$$
#'
#' Where the probability $\pi$ follows a $Beta$ distributions as:
#'
#' $$ \pi_{i,j,k} \sim Beta(\mu_{j,k},\phi)$$
#'
#' and
#'
#' $$ logit(\mu_{j,k}) = \alpha + \beta_{day}day_j + \beta_{assay}assay_k$$
#'
#' and
#'
#' $$ log(\phi_{j.k}) = \gamma + \delta_{day}day_j + \delta_{assay}assay_k
#'
#'
#' where $day_i$ has values of {0,1} for days 1 and 2, and $assay_i$ has
#' values of {0,1} for assays 1 and 2 respectively.
#'
#' The Odds Ratio (OR) for the probability of the day 2 vs day 1 is $e^{\beta_{day}}$
#'
#' the OR for the probability of assay 1 vs 2 is $e^{\beta_{assay}}$
#'
#' The relative change in dispersion parameter day 2 vs day 1 is $e^{\delta_{day}}$
#'
#' The relative change in dispersion parameter of  assay 2 vs 1 is $e^{\delta_{assay}}$
#'
#' The analysis is made using `rjags` which use shape parameters for the $Beta$
#' distribution. The conversion between precision and shape parameters is as follow:
#'
#' $$\mu = shape_1\(shape_1 +shape_2)$$
#' $$\phi = 1\(1 + shape_1 + shape_2) $$
#'
#' and
#'
#' $$ shape_1 = \mu(1-\phi)/\phi$$
#' $$ shape_2 = (1-\mu)(1-\phi)/\phi$$
#'
#' Non informative priors are given for the parameters
#' $$ \alpha \sim Normal(0,001)$$
#' $$ \beta_{day}\sim Normal(0,0.001)$$
#' $$ \beta_{assay} \sim Normal(0,0.001)$$
#' $$ \gamma \sim Normal(0,001)$$
#' $$ \delta_{day} \sim Normal(0,001)$$
#' $$ \delta_{assay} \sim Normal(0,001)$$
#'
#' The marginal correlation of the probabilities is estimated using the Pearson
#' correlation coefficient and reported as the mean value of the four individual
#' estimated correlations.
#'
#' The model is fit with 10000 simulations of the posterior in 4 independent
#' chains. Library `ggmcmc` is use to check the convergence of the model.
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

#'
#' ## Analysis
#'
#+ echo = T, message = F
# Simulate the data


# definition of the model
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
  alpha ~ dnorm(0, 0.001)
  bday ~ dnorm(0, 0.001)
  bassay ~ dnorm(0, 0.001)

  mu00 <- ilogit(alpha)
  mu01 <- ilogit(alpha + bassay)
  mu10 <- ilogit(alpha + bday)
  mu11 <- ilogit(alpha + bday + bassay)

  # Parameters that modify the dispersion in a log linear model
  gamma ~ dnorm(0, 0.001)
  dday ~ dnorm(0, 0.001)
  dassay ~ dnorm(0, 0.001)

  phi00 <- exp(gamma)
  phi10 <- exp(gamma + dday)
  phi01 <- exp(gamma + dassay)
  phi11 <- exp(gamma + dday + dassay)

  # Conversion from mu and phi to shape_1 and shape_2
  shape1_00 <- mu00*phi00
  shape1_01 <- mu01*phi01
  shape1_10 <- mu10*phi10
  shape1_11 <- mu11*phi11

  shape2_00 <- (1-mu00)*phi00
  shape2_01 <- (1-mu01)*phi01
  shape2_10 <- (1-mu10)*phi10
  shape2_11 <- (1-mu11)*phi11

  for (i in 1:Nsubj){

    # Probabilities for the binomial model
    prob_00[i] ~ dbeta(shape1_00, shape2_00)
    prob_01[i] ~ dbeta(shape1_01, shape2_01)
    prob_10[i] ~ dbeta(shape1_10, shape2_10)
    prob_11[i] ~ dbeta(shape1_11, shape2_11)

    # Binomial model
    pos_00[i] ~ dbinom(prob_00[i], mosq_00[i])
    pos_10[i] ~ dbinom(prob_10[i], mosq_10[i])
    pos_01[i] ~ dbinom(prob_01[i], mosq_01[i])
    pos_11[i] ~ dbinom(prob_11[i], mosq_11[i])

    # for correlation of probabilities
    standscore_assay_00_01[i] <- ((prob_00[i]-mean_00)/sd_00)*((prob_10[i]-mean_10)/sd_10)
    standscore_assay_10_11[i] <- ((prob_10[i]-mean_10)/sd_10)*((prob_11[i]-mean_11)/sd_11)
    standscore_assay_00_10[i] <- ((prob_00[i]-mean_00)/sd_00)*((prob_10[i]-mean_10)/sd_10)
    standscore_assay_01_11[i] <- ((prob_01[i]-mean_01)/sd_01)*((prob_11[i]-mean_11)/sd_11)

  }

  # Estimation of Pearson correlations based on the estimated probabilities
  # sd value is change if 0 when model is initialized
  mean_00 <- mean(prob_00[])
  mean_01 <- mean(prob_01[])
  mean_10 <- mean(prob_10[])
  mean_11 <- mean(prob_11[])
  sd_00 <- ifelse(sd(prob_00[]) == 0, 1,sd(prob_00[]))
  sd_01 <- ifelse(sd(prob_01[]) == 0, 1,sd(prob_01[]))
  sd_10 <- ifelse(sd(prob_10[]) == 0, 1,sd(prob_10[]))
  sd_11 <- ifelse(sd(prob_11[]) == 0, 1,sd(prob_11[]))
  corr_00_01 <- sum(standscore_assay_00_01[])/(Nsubj - 1)
  corr_10_11 <- sum(standscore_assay_10_11[])/(Nsubj - 1)
  corr_00_10 <- sum(standscore_assay_00_10[])/(Nsubj - 1)
  corr_01_11 <- sum(standscore_assay_01_11[])/(Nsubj - 1)
}
"
# Data example
set.seed(123456)

data_wider <- readxl::read_excel("_data/cvia_085_test.xlsx") %>%
  mutate(subjectid = seq(1,n())) %>%
  mutate(pos_00 = dsfa_bas_pos) %>%
  mutate(mosq_00 = dsfa_bas_pos + dsfa_bas_neg) %>%
  mutate(pos_01 = dmfa_bas_pos) %>%
  mutate(mosq_01 = dmfa_bas_pos + dmfa_bas_neg) %>%
  mutate(pos_10 = dsfa_fin_pos) %>%
  mutate(mosq_10 = dsfa_fin_pos + dsfa_fin_neg) %>%
  mutate(pos_11 = dmfa_fin_pos) %>%
  mutate(mosq_11 = dmfa_fin_pos + dmfa_fin_neg) %>%
  rename(subject = subjectid) %>%
  select(subject, starts_with("mosq"), starts_with("pos") )

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

# Simulate the posterior
latent_sim <-
  coda.samples(
    model = betabinom_jags,
    variable.names = c("alpha","bday","bassay","gamma","dday","dassay","corr_00_01","corr_10_11","corr_00_10","corr_01_11"),
    n.iter =  10000000,
    thin = 10)

#+ echo = T
# Extract and Analyze the posterior
latent_sum <-  summary(latent_sim)
eff_size <-  effectiveSize(latent_sim)
latent_chains <- ggs(latent_sim) %>%
  filter(Parameter %in% c("alpha","bday","bassay","gamma","dday","dassay"))

#' ## Results
res <- cbind(latent_sum[[1]], latent_sum[[2]], eff_size)
res[c("alpha","bassay","bday","gamma","dassay","dday","corr_00_01","corr_10_11","corr_00_10", "corr_01_11"), c("eff_size", "Mean", "SD","50%","2.5%","97.5%")]

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
