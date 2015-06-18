function result = ModAv(ensemb)
% Multi-model averaging function. Supported methods:
% EWA    : Equal Weights Averaging
% BGA    : Bates-Granger Averaging
% AICA   : Akaike's Information Criterion Averaging     
% BICA   : Bayes Information Criterion Averaging    
% GRA    : Granger-Ramanathan Averaging    
% BMA    : Bayesian Model Averaging
% MMA    : Mallows Model Averaging
% MMAd   : Mallows Model Averaging (lie on simplex)
% 
% Input fields
% Xcali  : m1 X n model prediction matrix for calibration period. m1 = time
% length of calibration period, n = model number.
% Ycali  : m1 X 1 observation matrix for calibration period.
% Xeval  : m2 X n model prediction matrix for evaluation period.
% Yeval  : m2 X 1 observation matrix for evaluation period.
% p      : 1 X n matrix of model parameter numbers.
% method : abbreviation of method.
%
% Output fields
% weight    : model weights.
% RMSEcali  : RMSE of calibration period. 
% RMSEeval  : RMSE of evaluation period.
% chain     : (BMA, MMA and MMAd) model weights in all iterations.
% obj       : (BMA, MMA and MMAd) objective function value for all iterations.
% sigma     : (BMA) optimized model sigma.
%
% Warning: Method BMA could be slow and unstable.

%% Unpack
[n,k] = size(ensemb.Xcali); % length of times series, number of models
Y = ensemb.Ycali; p = ensemb.p; method = ensemb.method;
Xeval = ensemb.Xeval; Yeval = ensemb.Yeval; t = size(Yeval,1);

%% Settings of DEMC
sigMaxBMA = sqrt(50); % Maximum sigma for BMA
sizeMMA = 50; lifeMMA = 800; mutMMA = 0.1;
sizeBMA = 100; lifeBMA = 20000; mutBMA = 0.1;

%% Bias correction
Xc = NaN(size(ensemb.Xcali));
result.bias = NaN(2,k);
for i=1:k
    x = [ensemb.Xcali(:,i),ones(n,1)];
    result.bias(:,i) = (x'*x)\x'*Y;
    Xc(:,i) = x*result.bias(:,i);
    Xeval(:,i) = [Xeval(:,i),ones(t,1)]*result.bias(:,i);
end
clear ensemb

%% Calculation of different methods
rmse = sqrt(mean((Xc-repmat(Y,[1,k])).^2)); % root mean square error
result.RMSEcali = rmse;
switch method
    %% Equal Weights Averaging
    case 'EWA' 
        result.weight = rmse.^0/numel(rmse);
    %% Bates-Granger Averaging    
    case 'BGA' 
        result.weight = 1./rmse.^2./sum(1./rmse.^2);
    %% Akaike's Information Criterion Averaging    
    case 'AICA' 
        I = n*log(rmse.^2) + n + 2*p;
        Iscale = I - min(I);
        result.weight = exp(-Iscale/2)/sum(exp(-Iscale/2));
    %% Bayes Information Criterion Averaging    
    case 'BICA' 
        I = n*log(rmse.^2) + n + p*log(n);
        Iscale = I - min(I);
        result.weight = exp(-Iscale/2)/sum(exp(-Iscale/2));
    %% Granger-Ramanathan Averaging    
    case 'GRA' 
        result.weight = ((Xc'*Xc)\Xc'*Y)';
    %% Optimazation-required methods   
    case {'BMA','MMA','MMAd'}
        pop.type = 'DEMC';
        pop.logFlag = true;
        %% Bayesian Model Averaging
        if strcmp(method,'BMA') 
            pop.size = sizeBMA;
            pop.life = lifeBMA;
            pop.mutation = mutBMA;
            pop.name = 'ModAv_BMA';
            f = @(sig) normpdf(repmat(Y,[1,k]),Xc,repmat(sig,[n,1]));
            pop.func = @(x)... % beta and sigma, half and half
                sum(log(sum(repmat(x(1:k),[n,1]).*f(x(k+1:end)),2)));
            pop.bound = [zeros(1,2*k);[ones(1,k),ones(1,k)*sigMaxBMA]];
            pop.simplex = [true(1,k),false(1,k)];
            DEMCout = DEMC(pop);
            result.sigma = DEMCout.best(k+1:end);
        %% Mallows Model Averaging
        else
            pop.size = sizeMMA;
            pop.life = lifeMMA;
            pop.mutation = mutMMA;
            pop.greedy = true;
            [~,idx] = max(p);
            S2 = rmse(idx)^2;
            pop.func = ...
                @(beta) -0.5*(sum((Y-Xc*beta').^2)+2*sum(beta.*p.*S2));
            %% Mallows Model Averaging (simplex version)
            if strcmp(method,'MMAd') % Weight that lie on the simplex
                pop.simplex = true(1,k);
                pop.bound = [zeros(1,k);ones(1,k)];
                pop.name = 'ModAv_MMAd';
            else
                pop.bound = [-ones(1,k);ones(1,k)];
                pop.name = 'ModAv_MMA';
            end
            DEMCout = DEMC(pop);
        end
        result.weight = DEMCout.best(1:k);
        result.chain = DEMCout.chain;
        result.obj = DEMCout.obj;
        
    otherwise
        warning('Method in function ModAv not recongized.');return
end

Ypred = sum(repmat(result.weight,[t,1]).*Xeval,2);
result.RMSEeval = sqrt(mean((Ypred-Yeval).^2));
end