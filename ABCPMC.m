function pop = ABCPMC(pop)
% Approximate Bayesian Computing-Population Monte Carlo method.
%
% Input fields
% obs         : 1 X d, d = dimension of target 
% bound       : 2 X n matrix, where n is the parameter dimension. Lower
% (first row) and upper(second row) bounds of parameter.
% funcPrior   : prior distribution function to generate parameter.
% funcSummar  : function to model target from parameter.
% funcDist    : function to calculate distance between modeled and observed
% target.
% epsl        : 1 X i epsilon matrix, where i is the number of generations.
% size        : population size
%
% Output the same strucuture with extra fields
% inds     : individuals (parameters) for all generations.
% summar   : modeled target at all generations.
% weight   : weight matrix for all generations.
% C        : Covariance matrix for all geneations.
% acRate   : overall acceptance rate.

%% Expand parameters
pop.life = numel(pop.epsl);
pop.dms = numel(pop.funcPrior());
pop.but = repmat(pop.bound(1,:),pop.size,1);
pop.top = repmat(pop.bound(2,:),pop.size,1);

%% Create empty matrixes
pop.inds = NaN(pop.size,pop.dms,pop.life);
pop.summar = NaN(pop.size,numel(pop.obs),pop.life);
pop.weight = NaN(pop.size,pop.life);
pop.C = NaN(pop.dms,pop.dms,pop.life); %Covariance matix
count = 0;

%% First iteration
for i = 1:pop.size
    flag = true;
    while flag
        theta = pop.funcPrior();
        Y = pop.funcSummar(theta);
        rho = pop.funcDist(Y);
        flag = rho > pop.epsl(1);
        count = count + 1;
    end
    pop.inds(i,:,1) = theta;
    pop.summar(i,:,1) = Y;
    pop.weight(i,1) = 1/pop.size;
    disp(['Generation 1, individual ' num2str(i)])
end
pop.C(:,:,1) = 2*cov(pop.inds(:,:,1));

%% Second and later iterations
q_d = NaN(pop.size,1); %probability of theta(j-1) to theta(j)
for j = 2:pop.life
    for i = 1:pop.size
        flag = true;
        while flag
            idx = randsample(1:pop.size,1,true,pop.weight(:,j-1));
            theta = pop.inds(idx,:,j-1);
            thetaP = mvnrnd(theta, pop.C(:,:,j-1));
            %% Check for out-off-bound
            if isfield(pop,'bound')
                thetaP = cage(thetaP,pop);
            end
            Y = pop.funcSummar(thetaP);
            rho = pop.funcDist(Y);
            flag = rho > pop.epsl(j);
            count = count + 1;
            disp(['Searching... ' num2str(count) '. ' num2str(((j-1)*...
                pop.size+i)/(pop.size*pop.life)*100) '% done.'])
        end
        pop.inds(i,:,j) = thetaP; 
        pop.summar(i,:,j) = Y;
        for u=1:pop.size 
            q_d(u) = mvnpdf(pop.inds(u,:,j-1),thetaP,pop.C(:,:,j-1));
        end
        bot = pop.weight(:,j-1)'*q_d;
        pop.weight(i,j) = 1./bot;
        disp(['Generation ' num2str(j) ', individual ' num2str(i)])
    end
    pop.C(:,:,j) = 2*cov(pop.inds(:,:,j));
end
pop.acRate = pop.size*pop.life/count;
