function result = MCMC(bead)
% Markov Chain Monte Carlo with Random-Walk Metropolis.
% 
% Input fields
% One of 'loc', 'funcPrior', 'bound'(in descending priority) must be
% specified. 
% life      : length of Markov Chain.
% func      : liklihood function. Calculates the liklihood from parameter.
% funcProp  : (optional) proposal distribution function. Make new proposal
% according to parameter. Default is normal distribution with sigma = 1. 
% loc       : (optional) initial parameter.
% funcPrior : (optional) prior function to generate initial parameter.
% bound     : (optional) 2 X n matrix, where n is the parameter dimension. 
% Lower(first row) and upper(second row) bounds of parameter.
% 
% Output fields
% chain     : m X n matrix, where m = length and n = parameter dimension.
% time      : m X 1 matrix, time stamp of each step.
%
% Warning   : Using 'bound' might break detailed balance and result in
% gibberish posterior. Not recommended.

%% Set default proposal distribution function to normal distribution.
if ~isfield(bead,'funcProp')
    bead.funcProp = @ (x) normrnd(x,1);
end

%% Retrive bead.loc and bead.dms (parameter dimension)
if ~isfield(bead,'loc')
    if isfield(bead,'funcPrior')
        bead.loc = bead.funcPrior();
    else
        bead.dms = size(bead.bound,2);
        bead.loc = rand(1,bead.dms)*(bead.bound(2,:)-bead.bound(1,:))+...
            bead.bound(1,:);
    end 
end
bead.dms = numel(bead.loc);

%% Create empty matrixes
rng('shuffle');
result.chain = NaN(bead.life, bead.dms);
result.time = NaN(bead.life, 1);

%% Initialize the first bead 
bead.p = bead.func(bead.loc);
bead.gen = 1;

%% Start Markov Chian and record
tic;
result.chain(bead.gen) = bead.loc;
result.time(bead.gen) = toc;
while bead.gen<=bead.life
    bead = RWM(bead);
    result.chain(bead.gen) = bead.loc;
    result.time(bead.gen) = toc;
end
end