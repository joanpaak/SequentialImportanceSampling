# If the support of the proposal distribution is too small compared to the
# posterior, the rejuvenation step can fail, or at least become computationally
# expensive. 
#
# This file demonstrates a situation in which that can happen.
#
# In the demonstration the option nu (which controles the heaviness of the 
# tails of the proposal distribution) is changed from small (heavy tails)
# to large (light tails). Pay attention to how the particle set gets "stuck"
# on around the first posterior when the tails are too light.
#
# The model used for this demo is a simple one-parameter Gaussian:
#
# y ~ normal(mu, 1.0)
# mu ~ normal(0, 10)
#
# We add two observations: 4 and -4.
#
# After the first observation the posterior is biased towards positive values.
# Because the prior is so wide, effective sample size gets small and the 
# the particle set is rejuvenated. If the tails of the proposal distribution
# are light, all of the proposals are generated near the first posterior, and
# because the posterior is relatively wide, acceptance rate of these proposals 
# is high. These factors lead to the posterior approximation being severely
# biased after this tep.
#
# When the negative observation is added to the model, the posterior should be
# at mu = 0, but because the posterior is "stuck" at the wrong place, the 
# estimate is too high.
#
# The problem can also be alleviated by
#   1) Doing importance sampling without rejuvenation. This does not address
#      particle degeneracy.
#   2) Increasing tempering/rejuvenation steps. This increases computational
#      load.
#

setwd("~/Desktop/SequentialImportanceSampling/")
source("sis_r6.R")

set.seed(1312)

draw_from_prior = function(n){
  return(cbind(rnorm(n, 0, 10)))
}

prior = function(theta){
  return(dnorm(theta[,1], 0, 10))
}

likelihood = function(y, theta){
  return(dnorm(y[1,1], theta[,1], 1.0))
}

nus = seq(50, 5, length.out = 10)
est_mu = rep(NaN, length(nus))

for(i in 1:length(nus)){
  sis = SIS$new(
    draw_from_prior = draw_from_prior,
    prior = prior,
    likelihood = likelihood,
    n_particles = 10000,
    opt = list(
      nu = nus[i],
      rejuvenation_limit = 0.75,
      tempering_limit = 0.50
    )
  )
  sis$add_observation(4)
  sis$add_observation(-4)
  
  est_mu[i] = sis$get_marginal_mus()
}

plot(nus, est_mu, type = "b", log = "x",
     bty = "l", ylab = "E[mu]", xlab = expression(nu))
abline(h = 0, lty = 3)

