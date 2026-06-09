# Student t distribution with fixed number of degrees of freedom and
# scale parameter.
#
# y ~ student_t(mu, s = 1, nu = 1)
# mu ~ student_t(mu = 0, s = 10, nu = 1)
#
# When n is small, the posterior will -- due to the sharpness of the 
# likelihood -- take on all kinds of funny shapes. When n = 2, the 
# posterior is often bi-modal.
# 

setwd("~/Desktop/SequentialImportanceSampling/")
source("sis_r6.R")

draw_from_prior = function(n){
  return(cbind(rt(n, 1) * 10))
}

prior = function(theta, n){
  return(dt(theta[,1] / 10, 1))
}

likelihood = function(y, theta){
  return(dt(y[1,1] - theta[,1], 1))
}

sis = SIS$new(
  draw_from_prior = draw_from_prior,
  prior = prior,
  likelihood = likelihood,
  n_particles = 1000
)

n = 3
y = cbind(rt(n, 1) + 2)

for(i in 1:n){
  sis$add_observation(y[i,])
}

theta = sis$get_iid_sample()
d_theta = 0.01 
grid_x = seq(min(theta), max(theta), d_theta)
grid_y = sis$posterior(cbind(grid_x), y)

hist(theta, prob = TRUE, breaks = 40)
points(grid_x, (grid_y / sum(grid_y)) / d_theta, type = "l")

