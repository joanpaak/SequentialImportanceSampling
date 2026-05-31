#
# A normal distribution model in which
#   y ~ normal(mu, sigma)
#

source("sis_r6.R")

gen_theta = cbind(c(4, 2))
N = 100
y = cbind(rnorm(N, gen_theta[1,], 1))

draw_from_prior = function(n){
  cbind(rnorm(n, 0, 10), rgamma(n, 3, 1))
}

prior = function(theta){
  dnorm(theta[,1], 0, 10) * dgamma(theta[,2], 3, 1)
}

likelihood = function(y, theta){
  return(dnorm(y[1,1], theta[,1], theta[,2]))
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

x11()
theta = sis$get_iid_sample()
plot(theta)

# Compare with grid approximation
delta_x = 0.01
grid_mu = seq(min(theta[,1]), max(theta[,1]), delta_x)
grid_sigma = seq(min(theta[,2]), max(theta[,2]), delta_x)
grid_y = sis$posterior(as.matrix(expand.grid(grid_mu, grid_sigma)), y)

grid_y = matrix(grid_y, 
                ncol = length(grid_mu), 
                nrow = length(grid_sigma),
                byrow = TRUE)

grid_y = grid_y / sum(grid_y) / delta_x

x11()
par(mfrow = c(1, 2))
hist(theta[,1], prob = TRUE)
points(grid_mu, colSums(grid_y) , type = "l")

hist(theta[,2], prob = TRUE)
points(grid_sigma, rowSums(grid_y) , type = "l")
