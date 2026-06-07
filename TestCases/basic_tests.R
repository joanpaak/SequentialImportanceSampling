#
# The SIS class is tested here by running it with various different 
# settings for a binomial model and then calculating the kl divergence
# to the true posterior. The divergences should, obviously, all be 
# reasonably close to zero.
#

setwd("~/Desktop/SequentialImportanceSampling/")
source("sis_r6.R")

# INPUT
#  theta : monte carlo sample
#  target : a function for calculating the target density
kl_divergence = function(theta, target){
  d = density(theta)
  dx = diff(d$x)[1]
  
  kl_div = sum(target(d$x) * log(target(d$x) / d$y) * dx, na.rm = TRUE)
  
  return(kl_div)
}

draw_from_prior = function(n){
  return(cbind(rbeta(n, 1, 1)))
}

prior = function(theta){
  return(dbeta(theta[,1], 1, 1))
}

likelihood = function(y, theta){
  dbinom(y[1,1], 1, theta[,1])
}

n = 40
y = cbind(rbinom(n, 1, 0.8))

sis_1 = SIS$new(
  draw_from_prior = draw_from_prior,
  prior = prior,
  likelihood = likelihood,
  n_particles = 5000,
  opt = list(
    logging = TRUE
  )
)

sis_2 = SIS$new(
  draw_from_prior = draw_from_prior,
  prior = prior,
  likelihood = likelihood,
  n_particles = 5000,
  opt = list(
    logging = TRUE,
    rejuvenation_limit = 1.2,
    tempering_limit = -1.2
  )
)

sis_3 = SIS$new(
  draw_from_prior = draw_from_prior,
  prior = prior,
  likelihood = likelihood,
  n_particles = 5000,
  opt = list(
    logging = TRUE,
    rejuvenation_limit = -1.2,
    tempering_limit = 1.2
  )
)

sis_4 = SIS$new(
  draw_from_prior = draw_from_prior,
  prior = prior,
  likelihood = likelihood,
  n_particles = 5000,
  opt = list(
    logging = TRUE,
    rejuvenation_limit = -1.2,
    tempering_limit = -1.2
  )
)

sis_5 = SIS$new(
  draw_from_prior = draw_from_prior,
  prior = prior,
  likelihood = likelihood,
  n_particles = 5000,
  opt = list(
    logging = TRUE,
    rejuvenation_limit = 1.2,
    tempering_limit = 1.2
  )
)

for(i in 1:n){
  sis_1$add_observation(y[i,,drop=FALSE])
  sis_2$add_observation(y[i,,drop=FALSE])
  sis_3$add_observation(y[i,,drop=FALSE])
  sis_4$add_observation(y[i,,drop=FALSE])
  sis_5$add_observation(y[i,,drop=FALSE])
}

target = function(x){
  dbeta(x, 1 + sum(y), 1 + n - sum(y))
}

cat(
  "These sould all be reasonably close to zero:",
  kl_divergence(sis_1$get_iid_sample(), target),
  kl_divergence(sis_2$get_iid_sample(), target),
  kl_divergence(sis_3$get_iid_sample(), target),
  kl_divergence(sis_4$get_iid_sample(), target),
  kl_divergence(sis_5$get_iid_sample(), target),
  sep = "\n"
)

