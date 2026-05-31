#
# A normal distribution model in which
#   y ~ normal(mu, 1.0)
#


source("sis_r6.R")

gen_theta = cbind(c(4))
N = 100
y = cbind(rnorm(N, gen_theta[,1], 1))

draw_from_prior = function(n){
  cbind(rnorm(n, 0, 10))
}

prior = function(theta){
  return(dnorm(theta[,1], 0, 10))
}

likelihood = function(y, theta){
  return(dnorm(y[1,1], theta[,1], 1.0))
}

sis = SIS$new(
  draw_from_prior = draw_from_prior,
  prior = prior,
  likelihood = likelihood,
  n_particles = 1e3
)

for(i in 1:N){
  sis$add_observation(y[i,,drop=FALSE])
}

theta = sis$get_iid_sample()
hist(theta, prob = TRUE)

# Compare with grid approximation
delta_x = 0.01
grid_x = seq(min(theta[,1]), max(theta[,1]), delta_x)
grid_y = sis$posterior(cbind(grid_x), y)
grid_y = grid_y / sum(grid_y) / delta_x

hist(sis$get_iid_sample(), prob = TRUE)
points(grid_x, grid_y , type = "l")

