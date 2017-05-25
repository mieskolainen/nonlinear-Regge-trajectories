% Mesonic cost function
% ------------------------------------------------------------------------
%
% Input:  x = d-dimensional parameter vector
%
% Output: y_w = weighted residual vector
%         y_p = pure residual vector
%        chi2 = chi2 cost
%
% mikael.mieskolainen@cern.ch

function [y_w, y_p, chi2] = Mcost(x, N_iters)

if (nargin == 1)
    N_iters = 15; % Default
end

% Parameters to be optimized
global sn;
global c;
global alpha0;
global lambda;
global init_lambda;

% Data
global N_M;
global N_W;
global N_M_err;
global N_W_err;
global N_J;
global N_J_err;

% Collect current parameter values
sn     = x(1:3);
c      = x(4:6);
alpha0 = x(7);

% First values based on linear trajectory
lambda = init_lambda(sn);

% First recursive fitting of lambda, based on current values
for k = 1:N_iters
    for n = 1:length(sn)
        % abs, real, due to numerical instability of recursion relation
        % lambda should be real and positive
        lambda(n) = abs(real(re_aM_s(sn(n))));
    end
end

% Model estimates
spin  = zeros(1,length(N_M));
width = zeros(1,length(N_M));
for i = 1:length(N_M)
    spin(i)  = re_aM_s(N_M(i)^2);
    width(i) = im_aM_s(N_M(i)^2) / (N_M(i) * re_apM_s(N_M(i)^2));
end

% Cost vector, because mass^2 approx~ spin,
% we use error on M^2 directly for the error on J
y_w = [(spin - N_J) ./ N_J_err, (width - N_W) ./ N_W_err];
y_p = [(spin - N_J), (width - N_W)];
chi2 = sum(y_w.^2);

end