# CATEGORICAL MODEL
#
# This example demonstrates a situation in which transformations
# are necessary for the rejuvenation step: it is almost certain that
# no proper unit simplexes are sampled. This, of course, requires
# the Jacobian of the transformation to be taken into account.
#
# Note that the functions for the transformations are not the most 
# efficient, focus is more on clarity.
#
# A Stan program is supplied for testing the convergence of the
# algorithm to the correct posterior.

source("sis_r6.R")

#### FUNCTIONS #####

# Dirichlet PDF and RNG from:
# https://en.wikipedia.org/wiki/Dirichlet_distribution
ddirichlet = function(x, alpha){
  (1.0 / (prod(gamma(alpha)) / gamma(sum(alpha)))) *
    prod(x^(alpha - 1))
}

rdirichlet = function(n, alpha){
  x = matrix(NaN, ncol = length(alpha), nrow = n)
  
  for(i in 1:n){
    x[i,] = rgamma(length(alpha), alpha, 1)
    x[i,] = x[i,] / sum(x[i,])
  }
  return(x)
}

# Transform vector into unit simplex.
# INPUT:
#   y  : vector on unconstrained scale
# OUPUT:
#  list with attributes:
#    theta : the transformed vector
#    J     : Jacobian of the transformation
transformToUnitVector = function(y){
  x = y
  K = length(x) + 1
  x[2:(K - 1)] = exp(x[2:(K - 1)])
  x = cumsum(x)
  
  p = rep(0, K)
  
  # The first and the last values:
  p[1] = plogis(x[1])
  p[K] = 1 - plogis(x[K-1])
  
  for(i in 2:(K - 1)){
    p[i] = plogis(x[i]) - plogis(x[i - 1])
  }
  
  J = prod(c(1, exp(y[2:(K-1)])) * (plogis(x) * (1 - plogis(x))))
  
  return(list(theta = p, J = J))
}

# Transformation from unit simplex to unconstrained
# space. Note that the transformation reduces the length
# of the vector by one,
#
# INPUT
#  x : a unit simplex
# OUTPUT
# a vector on the unconstrained scale.
inverse_transform = function(x){
  k = length(x) - 1
  x_ = qlogis(cumsum(x))[1:k]
  
  return(c(x_[1], log(x_[2:k] - x_[1:(k-1)])))
}

draw_from_prior = function(n){
  x = rdirichlet(n, c(3, 4, 5))
  theta = matrix(NaN, ncol = 2, nrow = n)
  
  for(i in 1:n){
    theta[i,] = inverse_transform(x[i,])
  }
  
  return(theta)
}

# Prior probability distribution
# Input is assumed to be on the unconstrained scale and
# it is then transformed into unit simplices and Jacobian
# adjustments are applied.
prior = function(theta_raw){
  theta = matrix(NaN, 
                 nrow = nrow(theta_raw), 
                 ncol = ncol(theta_raw) + 1)
  J = rep(NaN, nrow(theta_raw))
  prior_prob = rep(NaN, nrow(theta_raw))
  
  for(i in 1:nrow(theta_raw)){
    x = transformToUnitVector(theta_raw[i,])
    theta[i,] = x$theta
    J[i] = x$J
    prior_prob[i] = ddirichlet(theta[i,], c(3, 4, 5)) * J[i]
  }
  
  return(prior_prob)
}

likelihood = function(y, theta_raw){
  theta = matrix(NaN, nrow = nrow(theta_raw), ncol = ncol(theta_raw) + 1)
  
  for(i in 1:nrow(theta_raw)){
    x = transformToUnitVector(theta_raw[i,])
    theta[i,] = x$theta
  }
  
  return(theta[,y[1,1]])
}

#### INSTANTIATE IMPORTANCE SAMPLING AND ADD DATA ####

sis = SIS$new(
  draw_from_prior = draw_from_prior,
  prior = prior,
  likelihood = likelihood,
  n_particles = 1e3,
  opt = list(
    logging = TRUE
  )
)

gen_theta = rdirichlet(1, c(3, 4, 5))
n = 50
y = sample(1:length(gen_theta), n, TRUE, gen_theta)

for(i in 1:n){
  sis$add_observation(y[i])
}

theta = sis$get_iid_sample()
theta = t(apply(theta, 1, 
                function(x) transformToUnitVector(x)$theta))

matplot(t(theta), pch = 1, col = rgb(0, 0, 0, 0.05),
        ylim = c(0, 1))
points(c(1, 2, 3), gen_theta, pch = 19, col = "red")

sis$log$n_accepted

#### COMPARISON WITH STAN ####

cat_mod = rstan::stan_model(
  "ExampleModels/CategoricalModel/cat_model.stan")
cat_mod_fit = rstan::sampling(
  cat_mod,
  data = list(
    N = n,
    K = length(gen_theta),
    y = y,
    prior = c(3, 4, 5)
  ), chains = 2
)

theta_stan = as.matrix(cat_mod_fit)

x11()
par(mfrow = c(2, 3))
hist(theta_stan[,1], prob = TRUE)
points(density(theta[,1]), type = "l")

hist(theta_stan[,2], prob = TRUE)
points(density(theta[,2]), type = "l")

hist(theta_stan[,3], prob = TRUE)
points(density(theta[,3]), type = "l")

plot(theta_stan[,1:2])
points(theta[,1:2], col = "red")

plot(theta_stan[,2:3])
points(theta[,2:3], col = "red")

plot(theta_stan[,c(1, 3)])
points(theta[,c(1, 3)], col = "red")
