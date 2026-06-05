#
# Estimating the two-dimensional Rosenbrock function
#

setwd("~/Desktop/SequentialImportanceSampling/")
source("sis_r6.R")

rosenbrock = function(x, y, a = 1, b  = 100){
  return((a - x)^2 + b * ((y - x^2)^2))
}

x = seq(-2, 2, 0.1)
y = seq(-2, 4, 0.1)
z = matrix(NaN, ncol = length(x), nrow = length(y))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    z[j,i] = rosenbrock(x[i], y[j])
  }
}

contour(y, x, z)

# To make the estimation sequential, the Rosenbrock function can be 
# squeezed by raising it to power of (1.0 / PSEUDO_N) and then 
# "observations" are added given times. 

PSEUDO_N = 100

draw_from_prior = function(n){
  return(cbind(
    runif(n, -2, 2),
    runif(n, -2, 4)
  ))
}

prior = function(theta){
  dunif(theta[,1], -2, 2) * dunif(theta[,2], -2, 4)
}

likelihood = function(y, theta){
  rosenbrock(theta[,1], theta[,2])^(1 / PSEUDO_N)
}


sis = SIS$new(
  draw_from_prior = draw_from_prior,
  prior = prior,
  likelihood = likelihood,
  n_particles = 10000
)

for(i in 1:PSEUDO_N){
  sis$add_observation(matrix(NaN, nrow = 1))
}

theta = sis$get_iid_sample()

hist(theta[,1], prob = TRUE)
points(x, colSums(z) / sum(colSums(z)) / 0.1, type = "l")

hist(theta[,2], prob = TRUE)
points(y, rowSums(z) / sum(rowSums(z)) / 0.1, type = "l")

contour(x, y, t(z))
points(theta, pch = 19, col = rgb(0, 0, 0, 0.01))
