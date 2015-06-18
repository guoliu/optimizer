Optimizer tools written in Matlab. Including MCMC, DEMC, Ensemble Kalman filter, Approximate Bayesian Computing-Population Monte Carlo, and modeling averaging methods. Most algorithm based on class CEE 290 taught by Prof. [Jasper A. Vrugt](http://faculty.sites.uci.edu/jasper/).

All input and output of functions are single structures with multiple fields. Same fields name were used across different functions and can be piped together for testing purpose. 


## Documentation

#### function DEMC
Differential Evolution / Differential Evolution Markov Chain. In the case of **'DE'**, this function also plot a figure for the objective function of all generations.
##### Input fields
- func       : evaluation (objective) function.
- bound      : 2 X n matrix, where n is the parameter dimension. Lower (first row) and upper(second row) bounds of parameter.
- size       : size of population.
- life       : maximum generation of population.
- funcPrior  :(optional) prior distribution function to generate first
 generation. Default is uninformative prior within bounds.
- type       :(optional) 'DE' or 'DEMC'. The former seek the global minimum, the later maps the posterior. Default is 'DEMC'.
- greedy     :(optional) logical scalar, faster convergence but less stable. Breaks detailed balance and should not be used with 'DEMC'.
- logFlag    :(optional) logical scalar, use log liklihood or not. Default is 'false'.
- mutation   :(optional) scalar, mutation rate.
- CR         :(optional) scalar, crossover rate.
- simplex    :(optional) 1 X n logical maxtrix. 'true' specifies the parameters that should be positive and sum up to 1.
name       :(optional) name of output figure. Default 'DEMCtest'.

##### Output fields
- chain      : population at each generation.
- obj        : evaluation (objective) function value at each generation.
- best       : best individual at last generation.



#### function MCMC
Markov Chain Monte Carlo with Random-Walk Metropolis.
##### Input fields
One of 'loc', 'funcPrior', 'bound'(in descending priority) must be specified. 
- life      : length of Markov Chain.
- func      : liklihood function. Calculates the liklihood from parameter.
- funcProp  : (optional) proposal distribution function. Make new proposal
- according to parameter. Default is normal distribution with sigma = 1. 
- loc       : (optional) initial parameter.
- funcPrior : (optional) prior function to generate initial parameter.
- bound     : (optional) 2 X n matrix, where n is the parameter dimension. 
- Lower(first row) and upper(second row) bounds of parameter.
 
##### Output fields
- chain     : m X n matrix, where m = length and n = parameter dimension.
- time      : m X 1 matrix, time stamp of each step.

##### Warning   : Using 'bound' might break detailed balance and result in gibberish posterior. Not recommended.


#### function EnKF
Ensemble Kalman filter. Updated forecast and analysis states using observation and ensemble Kalman filter.

##### Input fields
- R        : measurement error (sigma).
- Fstates  : forecast states. 1 X n matrix, n = number of state
- Astates  : analysis states. 1 X n matrix.
- obs      : observation states. Scaler.
- func     : model operator, from current state to next state. 

##### Output 
The same structue with updated Fstates and Astates, as well as estimation of model error (C).



#### function ModAv
Multi-model averaging function. Supported methods:
- EWA    : Equal Weights Averaging
- BGA    : Bates-Granger Averaging
- AICA   : Akaike's Information Criterion Averaging     
- BICA   : Bayes Information Criterion Averaging    
- GRA    : Granger-Ramanathan Averaging    
- BMA    : Bayesian Model Averaging
- MMA    : Mallows Model Averaging
- MMAd   : Mallows Model Averaging (lie on simplex)
 
##### Input fields
- Xcali  : m1 X n model prediction matrix for calibration period. m1 = time
- length of calibration period, n = model number.
- Ycali  : m1 X 1 observation matrix for calibration period.
- Xeval  : m2 X n model prediction matrix for evaluation period.
- Yeval  : m2 X 1 observation matrix for evaluation period.
- p      : 1 X n matrix of model parameter numbers.
- method : abbreviation of method.

##### Output fields
- weight    : model weights.
- RMSEcali  : RMSE of calibration period. 
- RMSEeval  : RMSE of evaluation period.
- chain     : (BMA, MMA and MMAd) model weights in all iterations.
- obj       : (BMA, MMA and MMAd) objective function value for all iterations.
- sigma     : (BMA) optimized model sigma.

##### Warning: Method BMA could be slow and unstable.



#### function ABCPMC
Approximate Bayesian Computing-Population Monte Carlo method.

##### Input fields
- obs         : 1 X d, d = dimension of target 
- bound       : 2 X n matrix, where n is the parameter dimension. Lower
- (first row) and upper(second row) bounds of parameter.
- funcPrior   : prior distribution function to generate parameter.
- funcSummar  : function to model target from parameter.
- funcDist    : function to calculate distance between modeled and observed
- target.
- epsl        : 1 X i epsilon matrix, where i is the number of generations.
- size        : population size

#### Output 
The same strucuture with extra fields.
- inds     : individuals (parameters) for all generations.
- summar   : modeled target at all generations.
- weight   : weight matrix for all generations.
- C        : Covariance matrix for all geneations.
- acRate   : overall acceptance rate.


## Minimal Working Example
For DEMC and MCMC:
```
%% Finding global minimum with DE

% Test function: 2d Ackley's function
pop.func = @(x) -20*exp(-.2*sqrt(.5*sum(x.^2)))...
    -exp(.5*sum(cos(2*pi*x)))+exp(1)+20;
pop.bound = [-5,-5;5,5];

% Setup the DE settings
pop.size = 50; pop.life = 200; pop.type = 'DE';
result = DEMC(pop); 
disp(['Global minimum from DE: ' num2str(result.best)])


%% Mapping posterior with MCMC and DEMC

% Test function: wide spread mixture of two normal distributions: 1/3 of
% N(-5,1) and 2/3 of N(5,1). We expect a mean of 1.67 and a std of 4.82
pop.func = @(x) 1/3*normpdf(x,-5,1)+2/3*normpdf(x,5,1);
pop.bound = [-100;100]+1.67;

% Using DEMC. It should be able to sample both peaks.
pop.type = 'DEMC'; pop.life = 1500;
result = DEMC(pop); burnt = result.chain(pop.life/2:end,:);
disp(['Mean estimated by DEMC: ' num2str(mean(burnt(:)))])
disp(['Std estimated by DEMC: ' num2str(std(burnt(:)))])

% Using MCMC. It should only be able to sample one peak (-5 or 5).
pop.life = 10000; pop.loc = 0;
result = MCMC(pop); burnt = result.chain(pop.life/2:end,:);
disp(['Mean estimated by MCMC: ' num2str(mean(burnt))])
disp(['Std estimated by MCMC: ' num2str(std(burnt))])
```
