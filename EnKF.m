function ensemb = EnKF(ensemb)
% Ensemble Kalman filter. Updated forecast and analysis states using
% observation and ensemble Kalman filter.
%
% Input fields
% R        : measurement error (sigma).
% Fstates  : forecast states. 1 X n matrix, n = number of state
% Astates  : analysis states. 1 X n matrix.
% obs      : observation states. Scaler.
% func     : model operator, from current state to next state. 
% 
% Output the same structue with updated Fstates and Astates, as well as
% estimation of model error (C).

%% Update forecast states
N = numel(ensemb.Fstates);
Acell = num2cell(ensemb.Astates);
ensemb.Fstates = cellfun(ensemb.func,Acell);

%% Estimate model error
E = mean(ensemb.Fstates);
A = ensemb.Fstates-E;
ensemb.C = A*A.'/(N-1);

%% Perturb observation with measurement eror
D = ones(1,N)*ensemb.obs + normrnd(0,ensemb.R,[1,N]);

%% Kalman gain
K = ensemb.C/(ensemb.C+ensemb.R);

%% Update analysis states
ensemb.Astates = ensemb.Fstates + K*(D-ensemb.Fstates);