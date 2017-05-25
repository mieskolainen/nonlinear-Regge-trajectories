% Perturbated fminsearch, to avoid local minimas
%
% mikael.mieskolainen@cern.ch

function [x1, final_cost] = pertfminsearch(costfunc, x0, N_random, N_recursion, sigma)

% Nelder-Mead global optimization
min_cost = 1e99;
x_best = x0;
optim = optimset('MaxFunEvals', 500, 'Display','off');

for r = 1:N_random
    fprintf('.');
    
    % Perturbations
    x1 = x0 + randn(size(x0))*sigma.*abs(x0);
    [x1, final_cost] = fminsearch(costfunc, x1, optim);
    if (final_cost < min_cost)
        x_best = x1; min_cost = final_cost;
    end
end
x1 = x_best;
fprintf('\n');

% Finish with the recursive optimization, makes a large difference often
for i = 1:N_recursion
    fprintf('.');
    [x1, final_cost] = fminsearch(costfunc, x1, optim);
end
fprintf('\n\n');

end