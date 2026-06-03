
# Sequential Importance Sampling with Tempering

## Quick start

First, instantiate a new SIS object. To do that, you need three functions: 

1) draw_from_prior(n) 

This function should be defined in such a way that it returns a matrix of draws from the prior, columns corresponding to parameter values. 

2) prior(theta) 

Should return prior probability for a matrix of theta values, ie. the prior probability for each row of theta. 

3) likelihood(y, theta) 

This function should calculate likelihood for a single observation for a matrix of draws from the current posterior distribution and return a vector of likelihoods, each corresponding to the likelihood given a corresponding row of parameter values from the matrix of theta.

Posterior calculation is done internally by automatically combining the functions 2 and 3 from above listing. However, this function exists as the public method \$posterior and is thus at your mercy if you want to fiddle with it. E.g. I use it in some example to construct a grid approximation of the posterior to compare it with the Monte Carlo approximation.

During rejuvenation proposals are sampled from unconstrained space. In the case of models with highly constrained priors one might want to apply transformations so that the acceptance probability would be, well, acceptable. These transformations, if they are non-linear, require Jacobian adjustments. The basic idea is demonstrated in the Binomial example. 

## Example models

### Normal model

This folder contains two models.

The first model with known variance is...

y \~ normal(mu, 1.0)

mu \~ normal(0, 10)

...and in the second model both mu and sigma are assumed to be unknown:

y \~ normal(mu, sigma)

mu \~ normal(0, 10)

sigma \~ gamma(3, 1)

### Linear regression

Standard linear regression with gaussian errors.

y \~ normal(a + bx, sigma)

a \~ normal(0, 10)

b \~ normal(0, 10)

sigma \~ gamma(3, 1)

### Binomial model

This example is crafted as an example of how to include Jacobian adjustments to your model. The model in itself is a simple binomial model in which...

y \~ Binomial(theta)

...but the catch is that theta is non-linearly transformed such that

plogis(theta) \~ Beta(3, 3)

In order for the algorithm to converge to the correct posterior, the Jacobian has to be included in the model.
