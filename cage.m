function pop = cage(pop)
% Boundary and simplex constraint.
%
% Input field
% inds    : m X n maxtrix. m is the number of parameter combinations, n is
% the dimension of each combination.
% bound   : 2 X n matrix, where n is the parameter dimension. Lower(first
% row) and upper(second row) bounds of parameter.
% method  : (optional) one of 'reflect' (reflect back the over shoot
% distance), 'fold' (fold back the over shoot distance from the other
% bound), 'bound' (set to bound). Default is 'reflect'.
% simplex : (optional) 1 X n logical maxtrix. 'true' specifies the
% parameters that should be positive and sum up to 1.
%
% Output the same structure with constrained 'inds'.

%% Simplex requirement
if isfield(pop,'simplex')
    pop.inds(pop.inds<0) = 0; % set to 0
    pop.inds(:,pop.simplex) = pop.inds(:,pop.simplex)...
        ./repmat(sum(pop.inds(:,pop.simplex),2),...
        [1,sum(pop.simplex)]);
end

%% Expand bound to top and buttom
if ~(isfield(pop,'but')&&isfield(pop,'top'))
    m = size(pop.inds,1);
    pop.but = repmat(pop.bound(1,:),m,1);
    pop.top = repmat(pop.bound(2,:),m,1);
end

%% Boundary handeling. Default 'reflect'.
low = pop.inds < pop.but;
upp = pop.inds > pop.top;
if any([any(low),any(upp)]) %
    if ~isfield(pop,'method')
        pop.method = 'reflect';
    end
    
    switch pop.method
        case 'reflect' %reflect back the distance
            pop.inds(low)= 2 * pop.but(low) - pop.inds(low);
            pop.inds(upp)= 2 * pop.top(upp) - pop.inds(upp);
            
        case 'bound' %set to bound
            pop.inds(low)= pop.but(low);
            pop.inds(upp)= pop.top(upp);
            
        case 'fold' %circle back from the other bound
            pop.inds(low) = pop.top(low) - (pop.but(low) - pop.inds(low));
            pop.inds(upp) = pop.but(upp) + (pop.inds(upp) - pop.top(upp));
            
        otherwise
            warning('Method in function cage not recognized.'); return
    end
    %% Recursively avoid still out-of-bounds values
    pop = cage(pop);
end
end