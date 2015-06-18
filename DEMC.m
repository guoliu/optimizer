function result = DEMC(pop)
% Differential Evolution / Differential Evolution Markov Chain. In the 
% case of 'DE', this function also plot a figure for the objective function
% of all generations.
%
% Input fields
% func       : evaluation (objective) function.
% bound      : 2 X n matrix, where n is the parameter dimension. Lower
% (first row) and upper(second row) bounds of parameter.
% size       : size of population.
% life       : maximum generation of population.
% funcPrior  :(optional) prior distribution function to generate first
% generation. Default is uninformative prior within bounds.
% type       :(optional) 'DE' or 'DEMC'. The former seek the global
% minimum, the later maps the posterior. Default is 'DEMC'.
% greedy     :(optional) logical scalar, faster convergence but less
% stable. Breaks detailed balance and should not be used with 'DEMC'.
% logFlag    :(optional) logical scalar, use log liklihood or not. Default
% is 'false'.
% mutation   :(optional) scalar, mutation rate.
% CR         :(optional) scalar, crossover rate.
% simplex    :(optional) 1 X n logical maxtrix. 'true' specifies the
% parameters that should be positive and sum up to 1.
% name       :(optional) name of output figure. Default 'DEMCtest'.
%
% Output fields
% chain      : population at each generation.
% obj        : evaluation (objective) function value at each generation.
% best       : best individual at last generation.

%% Set defaul parameters
if ~isfield(pop,'type')
    pop.type = 'DEMC';
end
if ~isfield(pop,'name')
    pop.name = [pop.type,'test'];
end
if ~isfield(pop,'logFlag')
    pop.logFlag = false; %logFlag default false
end

%% Initialize parameters
pop.dms = size(pop.bound,2); %Parameter dimension
pop.but = repmat(pop.bound(1,:),pop.size,1);
pop.top = repmat(pop.bound(2,:),pop.size,1);
pop = rmfield(pop, 'bound');

%% Initialize population
pop.gen = 1; rng('shuffle');
if isfield(pop,'funcPrior')
    pop.inds = pop.funcPrior();
else
    pop.inds = rand(pop.size, pop.dms).*(pop.top-pop.but)+pop.but;
end
pop = cage(pop);

%% Evaluate first generation
pop.fitness = multiEval(pop.func, pop.inds);

%% Evolve population
result.chain = NaN(pop.life,pop.size,pop.dms);
result.obj = NaN(pop.life,pop.size);
while pop.gen<=pop.life
    result.chain(pop.gen,:,:) = pop.inds;
    result.obj(pop.gen,:) = pop.fitness;
    pop = evolve(pop);
end

%% Find best individual
fit = multiEval(pop.func, reshape(result.chain(end,:,:),[pop.size,pop.dms]));
if strcmp(pop.type,'DE')
    [~,idx] = min(fit);
    %% Plot the objective function
    figure;
    plot(1:pop.life,result.obj,'.');
    xlabel('Step'); ylabel('Objective Function');
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    saveas(gcf, [pop.name,'.png']);
    close(gcf);
else
    [~,idx] = max(fit);
end
result.best = reshape(result.chain(end,idx,:),[1,pop.dms]);
end