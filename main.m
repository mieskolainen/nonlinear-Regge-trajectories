% Complex valued non-linear Regge trajectories, based on:
%
%
% [1] R. Fiore, L. Jenkovszky, V. Magas, F. Paccanoni, A. Papa,
%     "Analytic Model of Regge Trajectories", Eur.Phys.J. A (2001)
%
% [2] R. Fiore, L. Jenkovszky, F. Paccanoni, A. Prokudin,
%     "Baryonic Regge trajectories with analyticity constraints",
%     Phys.Rev.D 70 (2004) 054003
%
%
% In general, input DATA should contain at least 4 hadrons on the trajectory.
% Otherwise #parameters >= #measurements. However, one can fit the
% subcritical case of 3 hadrons only too.
%
%
% mikael.mieskolainen@cern.ch, 2017

clear;
close all;

addpath('./src/'); % Our functions

% Number of fminsearch recursion trials and random perturbations
N_random = 20;
N_recursion = 10;
sigma = 0.04;

% Loop over particle families
tic;
diary('trajectory_output.txt');

for family = [1 2 3]

% Data from PDG
global N_M;
global N_W;
global N_M_err;
global N_W_err;
global N_J;
global N_J_err;

% Model parameters to be fitted
global sn;     % real 3-vector
global c;      % real 3-vector
global alpha0; % real
global delta;  % real

% Derived/calculated
global lambda;
global init_lambda;

% Different fermionic and bosonic families
if (family == 1)
    
    % DATA INPUT
    type  = 'fermionic';
    traj  = 'N';
    names = {'N(939)', 'N(1680)', 'N(2220)', 'N(2700)'};
    
    N_M     = [0.939, 1.680, 2.220, 2.700];  % Masses (GeV)
    N_W     = [1e-5, 0.130, 0.400, 0.350];   % Widths (GeV)
    N_M_err = [1e-6, 0.010, 0.050, 0.100];   % Mass errors (GeV)
    N_W_err = [1e-6, 0.010, 0.055, 0.050];   % Width errors (GeV)
    N_J     = [1/2, 5/2, 9/2, 13/2];         % Total angular momentum J
    
    % Model parameter initial values
    sn      = [1.2, 2.5, 12];
    c       = [0.5, 4.0, 4200];
    alpha0  = -0.40;
    delta   = -0.45;
    
elseif (family == 2)
    
    % DATA INPUT
    type  = 'bosonic';
    traj  = '\rho+a_2';
    names = {'rho(770)', 'a2(1320)', 'rho3(1690)', 'a4(2040)', 'rho5(2350)', 'a6(2450)'};
    
    N_M     = [0.77526, 1.3183, 1.6888, 1.995, 2.330, 2.450];
    N_M_err = [0.00025, 0.0005, 0.0021, 0.010, 0.035, 0.130];
    N_W     = [0.1478,  0.1050, 0.1610, 0.257, 0.400, 0.400];
    N_W_err = [0.0009,  0.0020, 0.0100, 0.025, 0.100, 0.250];
    N_J     = [1, 2, 3, 4, 5, 6];

    % Model parameter initial values
    sn      = [-0.4, 1.6, 14];
    c       = [0.09, 0.35, 18];
    alpha0  = 0.49;

elseif (family == 3)
    
    % DATA INPUT
    type  = 'bosonic';
    traj  = 'f';
    names = {'f2(1270)','f4(2050)','f6(2150)'};
    
    N_M     = [1.2755, 2.018, 2.469];
    N_M_err = [0.0010, 0.011, 0.029];
    N_W     = [0.1867, 0.237, 0.283];
    N_W_err = [0.0025, 0.018, 0.040];
    N_J     = [2, 4, 6];
    
    % Model parameter initial values
    sn      = [0.03, 1.4, 3.9];
    c       = [0.15, 0.16, 5];
    alpha0  = 0.8;
end

% Linear trajectory estimate on (M^2,J) plane
Y = N_J(:);
X = [ones(length(N_M),1) N_M(:).^2];
ALPHA = X \ Y; % Least Squares Linear regression without errors

% alpha_n = Re alpha(s_n), this is used in the fit as the initial
% estimate
init_lambda = @(x) ALPHA(1) + real(x)*ALPHA(2); % real() for fit safety
lambda = init_lambda(sn);

% "Effective" J error from the linear trajectory J = a0 + a'*M^2
% by error propagation (Taylor expansion). Thus we do not need to invert the
% non-linear trajectory which gives us J as a function of M^2 (not M^2 as a
% function of J).
N_J_err = sqrt( (2*ALPHA(2)*N_M).^2.*N_M_err.^2 );


%% FIT

fprintf('============================================================\n');
fprintf('Fitting the %s family: %s \n\n', type, traj);
fprintf('Linear (affine) parameters: a0 = %0.3f, a'' = %0.3f \n', ALPHA(1), ALPHA(2));

if (strcmp(type, 'fermionic'))

    % --------------------------------------------------------------------
    % Optimize parameters
    % First encapsulate parameters into vector x. Lambda iterated inside cost
    % function
    x0 = [sn, c, alpha0, delta];

    % Chi2 cost function
    Bsumcost = @(x) sum(Bcost(x).^2);
    [x1, cost] = pertfminsearch(Bsumcost, x0, N_random, N_recursion, sigma);
    
    % Update parameters by final call, IMPORTANT due to lambda!
    Bcost(x1);

    % Then non-linear least squares (for error estimation with output Jacobian)
    optim = optimset('MaxFunEvals', 500, 'Display','off','Algorithm','levenberg-marquardt');
    [x,~,~,~,~,~,J] = lsqnonlin(@Bcost, x1, [], [], optim);

    % Update parameters by final call, IMPORTANT due to lambda!
    [~,y_residual,chi2] = Bcost(x);
    
elseif (strcmp(type, 'bosonic'))
    
    % First encapsulate parameters into vector x. Lambda iterated inside cost
    % function
    x0 = [sn, c, alpha0];
    
    % Chi2 cost function
    Msumcost = @(x) sum(Mcost(x).^2);
    [x1, cost] = pertfminsearch(Msumcost, x0, N_random, N_recursion, sigma);
    
    % Update parameters by final call, IMPORTANT due to lambda!
    Mcost(x1);

    % Then non-linear least squares (for error estimation with output Jacobian)
    optim = optimset('MaxFunEvals', 500, 'Display','off','Algorithm','levenberg-marquardt');
    [x,~,~,~,~,~,J] = lsqnonlin(@Mcost, x1, [], [], optim);
    
    % Update parameters by final call, IMPORTANT due to lambda!
    [~,y_residual,chi2] = Mcost(x);
end



%% UNCERTAINTY ESTIMATION

% Linearized, asymptotic standard parameter errors
% see e.g., http://people.duke.edu/~hpgavin/ce281/lm.pdf

% Then normalize by NDF (#measurements - #parameters + 1), or by 1 if we
% have less measurements than parameters
ndf = size(J,1) - size(J,2) + 1;
w2 = y_residual(:)'*y_residual(:) / max(1, ndf);

fprintf('Fit results: \n\n');
fprintf('Chi2 / NDF = %0.3f / %d \n', chi2, ndf);

% A.) Construct error weight matrix based on empirical residuals
W = diag(1./w2);

% B.) Construct error weight matrix based on data errors
%err2 = [N_J_err, N_W_err].^2;
%W = diag(1./err2);

% Calculate covariance matrix for the parameters,
% pseudoinverse in a case we have undertermined Jacobian
% (= less measurements than parameters, or collinear/highly correlated measurements)
covar = pinv(full(J'*W*J));
x_err = sqrt(diag(covar));

% Collect parameter errors
sn_err = x_err(1:3);
c_err  = x_err(4:6);
alpha0_err = x_err(7);
if (strcmp(type,'fermionic'))
    delta_err  = x_err(8);
end

fprintf('\nParameter correlation matrix for %s -family: \n', traj);
fprintf('-------------------------------------------------------------\n');
corrcov(covar,1)

fprintf('\nObtained parameters for %s -family: \n', traj);
fprintf('-------------------------------------------------------------\n');
for n = 1:length(sn)
   fprintf('s(%d)   = %0.4f +- %0.4f \n', n, sn(n), sn_err(n));
end
for n = 1:length(c)
   fprintf('c(%d)   = %0.4f +- %0.4f \n', n, c(n),  c_err(n));
end
   fprintf('alpha0 = %0.4f +- %0.4f \n', alpha0, alpha0_err);
if (strcmp(type,'fermionic'))
   fprintf('delta  = %0.4f +- %0.4f \n', delta,  delta_err);
end

fprintf('\nCalculated: \n');
for n = 1:length(lambda)
   fprintf('lambda[%d] = %0.4f \n', n, lambda(n));
end

fprintf('\n');


%% Trajectory plot

figure;
s = -3:1e-2:9; % GeV^2
plottrajectory(s, N_M, N_M_err, N_J, alpha0, names, type, traj);
cmd = sprintf('print -dpdf ./figs/%s_trajectory.pdf', traj);
eval(cmd);


%% Width plot

figure;
M2 = linspace(0,8,1000); % GeV^2
plotwidth(M2, N_M, N_M_err, N_W, N_W_err, names, type, traj);
cmd = sprintf('print -dpdf ./figs/%s_width.pdf', traj);
eval(cmd);


%% Complex plane plot

figure;
s = -1:1e-2:9; % GeV^2
plotcomplexplane(s, type, traj);
cmd = sprintf('print -dpdf ./figs/%s_complexplane.pdf', traj);
eval(cmd);


end
diary off;
toc;