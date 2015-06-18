function pop = evolve(pop)
rng('shuffle');
grePcn = 90; % percentile of paerent in greedy method
lambda = 0.4; % strength of greedy jump
%% Define gamma
if mod(pop.gen,10)==0
    gamma = 1;
else
    gamma = 2.4/sqrt(2*pop.dms);
end

%% Generate parent matrix p, without repeating
s = (1:pop.size)';
p1 = mod(randi([0 pop.size-2],pop.size,1)+s,pop.size)+1;
p2 = randi([1,pop.size-2],pop.size,1);
p2 = (p2>=s|p2>=p1)+p2; p2 = (p2>=s&p2>=p1)+p2;

%% Greedy or gentle
if isfield(pop,'greedy') && (mod(pop.gen,10)~=0)
    starPool = find(pop.fitness>=prctile(pop.fitness,grePcn));
    ps = starPool(randi(numel(starPool),pop.size,1));
    pad = lambda*(pop.inds(ps,:)-pop.inds);
else
    pad = normrnd(0,1e-6,[pop.size,pop.dms]);
end

%% Create offspring embryo
embryo = pop.inds + gamma*(pop.inds(p1,:) - pop.inds(p2,:)) + pad;

%% Add crossover to offspring
if isfield(pop,'CR')
    trandition = rand(pop.size,pop.dms) > pop.CR;
    embryo(trandition) = pop.inds(trandition);
end

%% Add random (yet within range) mutations to offspring
if isfield(pop,'mutation')    
    weirdos = rand(pop.size,pop.dms) < pop.mutation;
    pool = rand(pop.size, pop.dms).*(pop.range(:,:,2)...
        - pop.range(:,:,1))+pop.range(:,:,1);
    embryo(weirdos) = pool(weirdos);
end

%% Reflect back outliers and meet simplex requirement
shell = pop; shell.inds = embryo;
shell = cage(shell); embryo = shell.inds;

%% Evaluate and substitue offspring
fitness = multiEval(pop.func, embryo);
switch pop.type
    case 'DEMC'
        if pop.logFlag
            snob = exp(fitness-pop.fitness)>rand(pop.size,1);
        else
            snob = fitness./pop.fitness>rand(pop.size,1);
        end
    case 'DE'
        snob = fitness<pop.fitness;
end
pop.inds(snob,:) = embryo(snob,:);
pop.fitness(snob,:) = fitness(snob,:);

%% Change generation number
pop.gen = pop.gen + 1;