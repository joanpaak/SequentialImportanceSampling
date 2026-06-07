#
# Test the do_rejuvenation() method.
#
# The idea here is to draw particles from a wide prior (uniform(-20, 20))
# while the posterior is norm(0, 1), then run the do_rejuvenation() method
# many times, while at the same time seeing if it converges to the target
# distribution. This is measured by calculating KL divergence between KDE
# of the particle set and the true density.
#
# If everything works, KL divergence should fall near zero and stay there.
#
# Note that the default posterior is overwritten and no  explicit 
# likelihood nor prior is used.
#

setwd("~/Desktop/SequentialImportanceSampling/")
source("sis_r6.R")

sis = SIS$new(
  draw_from_prior = function(n){
    cbind(runif(n, -20, 20))
  },
  
  prior = function(theta){
    return(0)
  },
  
  likelihood = function(y, theta){
    return(0)
  },
  
  n_particles = 1000
)

sis$posterior = function(theta, y){
  dnorm(theta[,1])
}

kl_divergence = c()

for(i in 1:20){
  sis$do_rejuvenation()
  theta = sis$get_iid_sample()
  d = density(theta)
  kl_divergence[i] = sum(dnorm(d$x) * log(dnorm(d$x) / d$y) * diff(d$x)[1])
}

plot(kl_divergence, type = "b",
     xlab = "Iteration", ylab = "KL divergence",
     bty = "l")

