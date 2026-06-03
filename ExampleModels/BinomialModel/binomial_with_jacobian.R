#
# Binomial model in which
#
#   y ~ Binomial(plogis(theta))
#   plogis(theta) ~ beta(3, 3)
#
# This demonstrates how to apply the Jacobian to take into account
# the non-linear transformation of theta.
# 
#
setwd("~/Desktop/Current Projects/SIS_REPO/")

source("sis_r6.R")

y = c(1, 1, 1)
n = length(y)

sis = SIS$new(
  draw_from_prior = function(n){
    cbind(qlogis(rbeta(n, 3, 3)))
  },
  prior = function(theta){
    dbeta(plogis(theta[,1]), 3, 3) * abs(dlogis(theta[,1]))
  },
  likelihood = function(y, theta){
    dbinom(y, 1, plogis(theta[,1]))
  },
  n_particles = 1e4
)


for(i in 1:n){
  sis$add_observation(y[i])
}

hist(plogis(sis$get_iid_sample()), prob = TRUE)
curve(dbeta(x, 3 + sum(y), 3 + (n - sum(y))), add = TRUE, col = "red")
