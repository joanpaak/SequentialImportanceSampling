# Standard linear regression with Gaussian errors 
#
# y[,2] ~ normal(theta[1] + theta[2] * y[,1], theta[3])
# theta[1] ~ normal(0, 10)
# theta[2] ~ normal(0, 10)
# theta[3] ~ gamma(3, 1)

source("sis_r6.R")

gen_theta = c(0, 1, 1)
N = 100
y = matrix(NaN, ncol = 2, nrow = N)
y[,1] = runif(N, -4, 4)
y[,2] = rnorm(N, gen_theta[1] + gen_theta[2] * y[,1], gen_theta[3])

draw_from_prior = function(n){
  cbind(rnorm(n, 0, 10), rnorm(n, 0, 10), rgamma(n, 3, 1))
}

prior = function(theta){
  dnorm(theta[,1], 0, 10) * dnorm(theta[,2], 0, 10) * 
    dgamma(theta[,3], 3, 1)
}

likelihood = function(y, theta){
  dnorm(y[,2], theta[,1] + theta[,2] * y[,1], theta[,3])
}

sis = SIS$new(draw_from_prior = draw_from_prior,
              prior = prior,
              likelihood = likelihood,
              n_particles = 2e3)

for(i in 1:N){
  sis$add_observation(y[i,])
}

theta = sis$get_iid_sample()

x11()
par(mfrow = c(1, 3))
hist(theta[,1]); abline(v = gen_theta[1])
hist(theta[,2]); abline(v = gen_theta[2])
hist(theta[,3]); abline(v = gen_theta[2])
