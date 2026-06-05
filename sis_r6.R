
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
      resampling_limit = 0.75,
      tempering_limit = 0.50,
      auto_adjust_k = TRUE,
      min_k = 1,
      max_k = 40,
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
      # NOTE: The proposal distribution might propose particles that
      #       are beyond the support of the prior distribution. In
      #       these cases useless warnings are generated which slow
      #       the script down. 
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
    
    add_observation = function(y, force_tempering = FALSE){
      if(is.null(y)){
        warning("y was NULL, assuming it is an empty observation")
        y = matrix(NaN, nrow = 1)
      }
      if(!is.matrix(y)) y = rbind(y)
      
      self$y <- rbind(self$y, y)
      p = self$likelihood(y, self$theta)
      w_updated = p * self$w
      w_updated = w_updated / sum(w_updated)
      
      n_eff = 1.0 / sum(w_updated^2)
      
      if(is.nan(n_eff)) n_eff = 0
      
      if(n_eff < (self$n_particles * self$opt$tempering_limit) |
         force_tempering){
        if(self$opt$auto_adjust_k){
          self$opt$k <- self$find_k(y) 
        }
        
        self$add_observation_with_tempering(y)
      } else if(n_eff < (self$n_particles * self$opt$resampling_limit)){
        self$w <- w_updated
        self$resample_and_move()
      } else {
        self$w <- w_updated
      }
      
      if(self$opt$logging){
        self$log$particle_set[[length(self$log$particle_set) + 1]] = 
          list(
            w = self$w,
            theta = self$theta
          )
        self$log$n_eff <- append(self$log$n_eff, 1.0 / sum(self$w^2))
      }
    },
    
    
    # Tries to find the value of k (number of tempering iterations)
    # that would lead to effective sample size not collapsing below
    # the option "tempering limit".
    find_k = function(y){
      if(is.null(ncol(y))) y = rbind(y)

      for(k in self$opt$min_k:self$opt$max_k){
        p = self$likelihood(y, self$theta)^(1.0 / k)
        w_updated = p * self$w
        w_updated = w_updated / sum(w_updated)
        
        n_eff = 1.0 / sum(w_updated^2)
        
        if(is.nan(n_eff)) return(self$opt$max_k)
        
        if(n_eff > (self$n_particles * self$opt$tempering_limit)){
          return(k)
        }
      }
      
      return(self$opt$max_k)
    },
    
    add_observation_with_tempering = function(y){
      if(is.null(ncol(y))) y = rbind(y)
      
      if(self$opt$logging) {
        self$log$k = rbind(self$log$k, 
                           c(nrow(self$y), self$opt$k))
      } 
      
      for(i in 1:self$opt$k){
        lh = self$likelihood(y, self$theta)^(1.0 / self$opt$k)
        self$w <- self$w * lh
        self$w <- self$w / sum(self$w)
        n_eff = (1.0 / sum(self$w^2))
        
        if(is.nan(n_eff)){
          self$resample_and_move()
          return()
        }
        
        if(n_eff < (self$n_particles * self$opt$resampling_limit)){
          self$resample_and_move()
        }
      }
    },
    
    resample_and_move = function(){
      # TODO: Implement stratified resampling?
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
        warning(paste("Rejuvenation failed, only", 
                      length(inds_accepted), "accepted proposals"))
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

