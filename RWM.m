function bead = RWM(bead)

%% Create proposal 
loc = bead.funcProp(bead.loc);
%% Check bound
if isfield(bead,'bound')
    pop.inds = loc; pop.bound = bead.bound;
    pop = cage(pop);
    loc = pop.inds;
end

%% Calculate and compare Metropolis ratio
p = bead.func(loc);
if p/bead.p > rand()
    bead.loc = loc;
    bead.p = p;
end

%% Register sequence
bead.gen = bead.gen+1;