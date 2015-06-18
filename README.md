Optimizer tools written in Matlab. Including MCMC, DEMC, Ensemble Kalman filter, Approximate Bayesian Computing-Population Monte Carlo, and modeling averaging methods. Most algorithm based on class CEE 290 taught by Prof. [Jasper A. Vrugt](http://faculty.sites.uci.edu/jasper/).

function DEMC
=============
Differential Evolution / Differential Evolution Markov Chain. In the case of 'DE', this function also plot a figure for the objective function of all generations.

Input fields
------------
func       : evaluation (objective) function.
bound      : 2 X n matrix, where n is the parameter dimension. Lower (first row) and upper(second row) bounds of parameter.
size       : size of population.
life       : maximum generation of population.
funcPrior  :(optional) prior distribution function to generate first
generation. Default is uninformative prior within bounds.
type       :(optional) 'DE' or 'DEMC'. The former seek the global minimum, the later maps the posterior. Default is 'DEMC'.
greedy     :(optional) logical scalar, faster convergence but less stable. Breaks detailed balance and should not be used with 'DEMC'.
logFlag    :(optional) logical scalar, use log liklihood or not. Default is 'false'.
mutation   :(optional) scalar, mutation rate.
CR         :(optional) scalar, crossover rate.
simplex    :(optional) 1 X n logical maxtrix. 'true' specifies the parameters that should be positive and sum up to 1.
name       :(optional) name of output figure. Default 'DEMCtest'.

Output fields
------------
chain      : population at each generation.
obj        : evaluation (objective) function value at each generation.
best       : best individual at last generation.
