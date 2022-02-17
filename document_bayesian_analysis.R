#' ---
#' title: Bayesian model for the CVIA-085 trial with both precision and Mean change
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
#' $$ \pi_{i,j,k} \sim Beta(\mu_{j,k},\phi_{j,k})$$
#' 
#' and
#' 
#' $$ logit(\mu_{j,k}) = \alpha + \beta_{day}day_j + \beta_{assay}assay_k$$
#' 
#' and
#' 
#' $$ log(\phi_{j,k}) = \gamma + \delta_{day}day_j + \delta_{assay}assay_k$$
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
#' $$\mu = shape_1/(shape_1 +shape_2)$$
#' $$\phi = shape_1 + shape_2 $$
#' 
#' and
#' 
#' $$ shape_1 = \mu\phi$$
#' $$ shape_2 = (1-\mu)\phi$$
#' 
#' The intraclass correlation is: $ICC=1/(\phi-1)
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
