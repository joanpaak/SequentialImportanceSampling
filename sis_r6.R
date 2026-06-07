# R6 class for sequential importance sampling with optional rejuvenation steps
# and tempering.
#
# CALLING THE CONSTRUCTOR
#
# MANDATORY INPUTS:
# draw_from_prior : function with one parameter n, should return a matrix of 
#                   draws from the prior with parameters on columns.
# prior           : function for calculating the prior probability density for 
#                   a matrix of draws from the prior. Should return a vector 
#                   containing p(theta) for each row of the supplied matrix of 
#                   draws.
# likelihood      : function for calculating likelihood of a single observation 
#                   given a matrix of draws from the posterior. Should return a 
#                   vector containing p(y|theta) for each row of the supplied 
#                   matrix of draws.
# n_particles     : integer, the number of particles to use.
#
# OPTIONAL INPUTS:
# opt             : List of options to set. The settable options are:
# opt$k           : integer, a pre-set number of tempering steps.
# opt$auto_adjust_k : boolean, should k be adjusted automatically.
# opt$min_k         : integer, minimum number of tempering steps if 
#                     auto-adjusting is used.
# opt$max_k         : integer, maximum number of tempering steps if 
#                     auto-adjusting is used.
# opt$rejuvenation_limit : number between 0 and 1. If n_eff/n_particles falls 
#                          below this, a rejuvenation step is performed.
# opt$tempering_limit  : number between 0 and 1. if n_eff/n_particles falls  
#                        below this, the observation is added with tempering.
# opt$logging          : boolean, should logging be done. 
# opt$n_rejuvenation_steps : How many times the rejuvenation step is performed. 
#
# METHODS
# 
# add_observation(y) : adds a single new observation. If your observations are vectors, 
#                      e.g. in linear regression observation pairs, each row of y
#                      should correspond to a single observation. Set the optional argument
#                      force_tempering to TRUE if you want to... force tempering!
#                      Set the optional argument force_rejucenation to... you guessed it.
#                      Adds the observation to an internal data matrix that's referenced 
#                      during rejuvenation.
# get_iid_sample()   : samples from the current posterior using the current weights
# get_marginal_mus() : returns marginal means calculated using current weights
# get_marginal_sds() : returns marginal sds calculated using current weights
#
# The following methods are typically used internally but can 
# be useful for debugging or otherwise:
# do_rejuvenation()   : resamples using multinomial resampling, generates 
#                       proposals and accepts them according to their posterior 
#                       probability.
# do_tempering(y) : Same as add_observation but with tempering. NOTE: 
#                   DOES NOT ADD THE OBSERVATION TO THE INTERNAL DATA MATRIX
#                   THAT IS USED FOR CALCULATING ACCEPTANCE PROBABILITIES.
#                   Use add_observation with the argument force_tempering = 
#                   TRUE if you want to make sure tempering is used.
# posterior(theta, y) : calculation of posterior probabilities. 
#
# EXAMPLE:
# Gaussian model with known variance...
#
# y ~ normal(mu, 1.0)
# mu ~ normal(0, 1)
#
# ...would be set up like this:
#
# sis = SIS$new(
#   draw_from_prior = function(n){
#     return(cbind(rnorm(n)))
#   },
#   prior = function(theta){
#     return(dnorm(theta[,1]))
#   },
#   likelihood = function(y, theta){
#     return(dnorm(y[1,1], theta[,1], 1.0))
#   },
#   n_particles = 1000
#)
#
# Then we can add observations...
#
# for(i in 1:10){
#   sis$add_observation(rnorm(1))
# }
#
# ...and plot the posterior:
#
# hist(sis$get_iid_sample())
#
# TODO:
# - Improve logging
# - Detect underflow when tempering with large k

SIS <- R6::R6Class(
  "SIS",
  public = list(
    y     = NULL,
    theta = NULL,
    n_dim = NULL,
    w     = NULL,
    n_particles = NULL,
    likelihood  = NULL,
    prior       = NULL,
    posterior   = NULL,
    
    opt = list(
      k = 5,
      rejuvenation_limit = 0.75,
      tempering_limit = 0.50,
      auto_adjust_k = TRUE,
      min_k = 2,
      max_k = 40,
      n_rejuvenation_steps = 1,
      logging = FALSE
    ),
    
    log = list(
      n_eff = c(),
      k = matrix(NaN, ncol = 2,  nrow = 0),
      particle_set = list(),
      n_accepted = matrix(NaN, ncol = 2,  nrow = 0)
    ),
    
    initialize = function(draw_from_prior,
                          prior,
                          likelihood,
                          n_particles,
                          opt = NULL){
      if(!is.null(opt)){
        n = names(opt)
        for(i in 1:length(opt)){
          self$opt[[n[i]]] = opt[[i]]
        }
      }
      
      self$n_particles <- n_particles
      self$theta <- draw_from_prior(self$n_particles)
      self$n_dim <- ncol(self$theta)
      self$likelihood <- likelihood
      self$prior <- prior
      
      self$posterior  <- function(theta, y){
        post_prob = self$prior(theta)
        non_zero_inds = which(post_prob > 0)
        
        if(!is.null(self$y)){
          for(i in 1:nrow(y)){
            post_prob[non_zero_inds] = 
              post_prob[non_zero_inds] * 
              self$likelihood(y[i,,drop=FALSE], theta[non_zero_inds,,drop=FALSE])
          }
        }
        
        return(post_prob)
      }
      
      self$w <- rep(1.0 / self$n_particles, self$n_particles)
      self$y <- c()
      colnames(self$log$k) = c("t", "k")
      colnames(self$log$n_accepted) = c("t", "N accepted")
    },
    
    # Calculates the updated weights given an observation. Mainly a 
    # convenience method for doing error correction/emitting warnings
    # in a single place.
    #
    # INPUT
    #  y : a row of a data matrix
    #  k : (optional) tempering constant for the likelihood
    # OUTPUT
    #  result$w_updated : the updated weights
    #  result$n_eff     : the proportional effective sample size
    get_updated_weights = function(y, k = 1){
      w_updated = self$w * likelihood(y, self$theta)^(1/k)
      w_updated = w_updated / sum(w_updated)
      
      if(any(is.nan(w_updated))){
        warning("NaN values when updating weights, ",
                "replacing them with zeroes.")
        w_updated[is.nan(w_updated)] = 0
      }
      
      n_eff = (1.0 / sum(w_updated^2)) / self$n_particles
      
      return(list(
        w_updated = w_updated,
        n_eff = n_eff
      ))
    },
    
    add_observation = function(y, 
                               force_rejuvenation = FALSE, 
                               force_tempering = FALSE){
      if(is.null(y)) stop("y was null")
      if(!is.matrix(y)) y = rbind(y)
      
      ## Add observation to the internal data matrix
      self$y <- rbind(self$y, y)
      
      ## Update particle set
      update_result = self$get_updated_weights(y)
      
      n_eff_above_rejuvenation = update_result$n_eff > self$opt$rejuvenation_limit
      n_eff_above_tempering = update_result$n_eff > self$opt$tempering_limit
      
      if(force_tempering) n_eff_above_tempering = FALSE
      if(force_rejuvenation) n_eff_above_rejuvenation = FALSE
      
      if(n_eff_above_rejuvenation & n_eff_above_tempering){
        self$w <- update_result$w_updated
      } else if(n_eff_above_rejuvenation & !n_eff_above_tempering){
        self$do_tempering(y)
      } else if(!n_eff_above_rejuvenation & n_eff_above_tempering){
        self$w <- update_result$w_updated
        for(j in 1:self$opt$n_rejuvenation_steps){
          self$do_rejuvenation()
        }
      } else if(!n_eff_above_rejuvenation & !n_eff_above_tempering){
        self$do_tempering(y, force_rejuvenation = TRUE)
      }
      
      if(self$opt$logging){
        self$log$particle_set[[nrow(self$y + 1)]] = list(
          theta = self$theta,
          w = self$w
        )
      }
    },
    
    # Tries to find the value of k (number of tempering iterations)
    # that would lead to effective sample size not collapsing below
    # the option "tempering limit".
    adjust_k = function(y){
      for(k in self$opt$min_k:self$opt$max_k){
        update_result = self$get_updated_weights(y, k)
        
        if(update_result$n_eff > self$opt$tempering_limit){
          return(k)
        }
      }
      
      warning("Maximum iterations reached for finding k. ",
              "Proportional n_eff after updating weights ", update_result$n_eff)
      return(self$opt$max_k)
    },
    
    # Updates the weights using tempering and rejuvenation steps
    # whenever n_eff falls below rejuvenation_limit during the 
    # tempering. NOTE: Does not add y to the internal data matrix!
    do_tempering = function(y, force_rejuvenation = FALSE){
      if(self$opt$auto_adjust_k) self$opt$k = self$adjust_k(y)
      if(self$opt$logging){
        self$log$k <- rbind(
          self$log$k,
          c(nrow(self$y), self$opt$k)
        )
      }
      
      lh = self$likelihood(y, self$theta)
      
      for(i in 1:self$opt$k){
        self$w <- self$w * lh^(1/self$opt$k)
        self$w = self$w / sum(self$w)
        n_eff = (1.0 / sum(self$w^2)) / self$n_particles
        
        if(n_eff < self$opt$rejuvenation_limit | force_rejuvenation){
          for(j in 1:self$opt$n_rejuvenation_steps){
            self$do_rejuvenation()
          }
          lh = self$likelihood(y, self$theta)^(1/self$opt$k)
        }
      }
    },
    
    # Performs the rejuvenation step. First resamples particles by using
    # multinomial sampling, then generates proposals from a multidimensional
    # Gaussian with means and standard deviations corresponding to the current
    # posterior and then accepts/rejects those proposals using and independent
    # Metropolis-Hastings kernel.
    # 
    # The internal data matrix is used when calculating acceptance proposals.
    do_rejuvenation = function(){
      # TODO: Implement stratified rejuvenation?
      rs_inds = sample(1:self$n_particles, self$n_particles, TRUE, self$w)
      theta_rs = self$theta[rs_inds,,drop = FALSE]
      theta_prop = matrix(NaN, ncol = self$n_dim, nrow = self$n_particles)
      
      marg_mus = rep(NaN, self$n_dim)
      marg_sds = rep(NaN, self$n_dim)
      
      for(i in 1:self$n_dim){
        marg_mus[i] = sum(self$theta[,i] * self$w)
        marg_sds[i] = sqrt(sum((self$theta[,i] - marg_mus[i])^2 * self$w))
        
        theta_prop[,i] = rnorm(self$n_particles, marg_mus[i], marg_sds[i])
      }
      
      p_curr = rep(NaN, self$n_particles)
      p_prop = rep(NaN, self$n_particles)
      mh_ratio = rep(NaN, self$n_particles)
      
      p_curr = self$posterior(theta_rs[,,drop=FALSE], self$y[,,drop=FALSE])
      p_prop = self$posterior(theta_prop[,,drop=FALSE], self$y[,,drop=FALSE])
      
      for(i in 1:self$n_particles){
        # NOTE TO SELF: Remember that this is not the MH-kernel but
        # _independent_ MH kernel since the proposals are generated
        # from an independent distribution
        mh_ratio[i] = prod(dnorm(theta_rs[i,], marg_mus, marg_sds)) / 
          prod(dnorm(theta_prop[i,], marg_mus, marg_sds))
      }
      
      p_accept = (p_prop / p_curr) * mh_ratio
      p_accept[is.na(p_accept)] = 0
      
      s = runif(self$n_particles, 0, 1)
      inds_accepted = which(s < p_accept)
      self$theta <- theta_rs
      self$theta[inds_accepted, ] <- theta_prop[inds_accepted,,drop = FALSE]
      self$w <- rep(1.0 / self$n_particles, self$n_particles)
      
      if(length(inds_accepted) < 10){
        warning("Rejuvenation failed, only", 
                 length(inds_accepted), "accepted proposals")
      }
      
      if(self$opt$logging){
        self$log$n_accepted = rbind(
          self$log$n_accepted, 
          c(nrow(self$y), length(inds_accepted)))
      }
    },
    
    get_marginal_mus = function(){
      mus = rep(NaN, self$n_dim)
      
      for(i in 1:self$n_dim){
        mus[i] = sum(self$theta[,i] * self$w)
      }
      
      return(mus)
    },
    
    get_marginal_sds = function(){
      mus = self$get_marginal_mus()
      sds = rep(NaN, self$n_dim)
      
      for(i in 1:self$n_dim){
        sds[i] = sqrt(sum((self$theta[,i] * self$w)^2))	
      }
      
      return(sds)
    },
    
    get_iid_sample = function(){
      rs_inds = sample(1:self$n_particles, self$n_particles, TRUE, self$w)
      theta_rs = self$theta[rs_inds,,drop = FALSE]
      
      return(theta_rs)
    }
  )
)
