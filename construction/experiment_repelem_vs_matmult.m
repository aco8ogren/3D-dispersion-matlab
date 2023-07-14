clear; close all;

N = 1000000;

N_trials = 100;

E = rand(N,1);
nu = rand(N,1);

% Matrix multiplication
disp('matrix multiplication')
tic
for i = 1:N_trials
[(32-48*nu).*E 6*E 6*E -8*E (24*nu-6).*E (24*nu-6).*E 4*E (-24*nu+6).*E 3*E (-10+12*nu).*E -6*E (12*nu-3).*E 4*E 3*E (-24*nu+6).*E (-10+12*nu).*E (12*nu-3).*E -6*E (-4+12*nu).*E (-12*nu+3).*E (-12*nu+3).*E (-8+12*nu).*E -3*E -3*E];
end
toc

% Repelem
disp('repelem')
tic
s = size(E,1);
for i = 1:N_trials
E.*[32-48*nu repelem(6,s,1) repelem(6,s,1) repelem(-8,s,1) 24*nu-6 24*nu-6 repelem(4,s,1) -24*nu+6 repelem(3,s,1) -10+12*nu repelem(-6,s,1) 12*nu-3 repelem(4,s,1) repelem(3,s,1) -24*nu+6 -10+12*nu 12*nu-3 repelem(-6,s,1) -4+12*nu -12*nu+3 -12*nu+3 -8+12*nu repelem(-3,s,1) repelem(-3,s,1)];
end
toc

% Multiply by ones
disp('multiply by ones')
tic
s = size(E,1);
for i = 1:N_trials
E.*[32-48*nu 6*ones(s,1) 6*ones(s,1) -8*ones(s,1) 24*nu-6 24*nu-6 4*ones(s,1) -24*nu+6 3*ones(s,1) -10+12*nu -6*ones(s,1) 12*nu-3 4*ones(s,1) 3*ones(s,1) -24*nu+6 -10+12*nu 12*nu-3 -6*ones(s,1) -4+12*nu -12*nu+3 -12*nu+3 -8+12*nu -3*ones(s,1) -3*ones(s,1)];
end
toc

% Multiply by predefined ones
disp('multiply by predefined ones')
tic
o = ones(s,1);
s = size(E,1);
for i = 1:N_trials
E.*[32-48*nu 6*o 6*o -8*o 24*nu-6 24*nu-6 4*o -24*nu+6 3*o -10+12*nu -6*o 12*nu-3 4*o 3*o -24*nu+6 -10+12*nu 12*nu-3 -6*o -4+12*nu -12*nu+3 -12*nu+3 -8+12*nu -3*o -3*o];
end
toc

