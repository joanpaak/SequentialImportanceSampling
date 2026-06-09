
# Sequential Importance Sampling with Tempering

An R6 class for performing sequential importance sampling for approximating Bayesian posterior distiributions with weighted random numbers (particles). 

In sequential importance sampling observations are added to the algorithm one-by-one and during each step the weights are updated using the likelihood. This will often lead to only a few particles having significant weight. This phenomenon, known as particle degeneracy, is countered by adding rejuvenation steps into the algorithm.

During these rejuvenation steps proposed new particles are generated from a multidimensional t-distribution with means and variances corresponding to the current posterior. Acceptance probabilities are calculated using an independent Metropolis-Hastings kernel. Given that the posterior is sufficiently good approximation of the true posterior, a single application of this process should be enough to rejuvenate the particles.

In some models, however, a single observation might collapse the weights in such a way that the rejuvenation step is unable to "reset" the particles. In these cases I use tempering. 

During tempering, the likelihood is raised to some fractinal power, L(y | theta)^(1/k) and applied k times. Rejuvenation steps are performed inside this tempering loop.

## Quick start

First, instantiate a new SIS object. To do that, you need three functions: 

1) draw_from_prior(n) 

This function should be defined in such a way that it returns a matrix of draws from the prior, columns corresponding to parameter values. 

2) prior(theta) 

Should return prior probability for a matrix of theta values, ie. the prior probability for each row of theta. 

3) likelihood(y, theta) 

This function should calculate likelihood for a single observation for a matrix of draws from the current posterior distribution and return a vector of likelihoods, each corresponding to the likelihood given a corresponding row of parameter values from the matrix of theta.

Posterior calculation is done internally by automatically combining the functions 2 and 3 from above listing. However, this function exists as the public method \$posterior and is thus at your mercy if you want to fiddle with it. E.g. I use it in some example to construct a grid approximation of the posterior to compare it with the Monte Carlo approximation.

During rejuvenation proposals are sampled from unconstrained space. In the case of models with highly constrained priors one might want to apply transformations so that the acceptance probability would be, well, acceptable. These transformations, if they are non-linear, require Jacobian adjustments. The basic idea is demonstrated in the Binomial example, and a more realistic example is given in the categorical model example; it is more realistic in the sense that here an approprioate transformation is crucial for the rejuvenation step to be succesful.

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

### Categorical model

This example demonstrates non-linear transformations and how to apply Jacobian adjustments. 

The model is...

y \~ categorical(phi)

phi \~ dirichlet(3, 4, 5)

Because during the rejuvenation step proposals are drawn from a multidimensional normal distribution, it is almost certain that none of the proposals are unit simplexes, and so the prior probability of each proposal would by default be zero. 


