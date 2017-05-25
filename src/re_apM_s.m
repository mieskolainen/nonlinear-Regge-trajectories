% Mesonic trajectory, real part of the trajectory slope, Re alpha'(s)
% ------------------------------------------------------------------------
%
%     "Analytic Model of Regge Trajectories"
%     R. Fiore, L. Jenkovszky, V. Magas, F. Paccanoni, A. Papa,
%     Eur.Phys.J. A (2001)
%
%
% mikael.mieskolainen@cern.ch

function out = re_apM_s(s)

global c;

% Total number of resonances
N = length(c);

out = 0;
for n = 1:N
   out = out + c(n)*A(s,n);
end
out = out / sqrt(pi);

% Take real part due to floating points
out = real(out);

end

% Amplitudes
function out = A(s,n)

global sn;
global lambda;

out = 0;

if (theta(sn(n) - s) > 0 && abs(s/sn(n)) < (1 + eps(s))) % Check sanity
out = out + (gamma(lambda(n)+3/2)) / (gamma(lambda(n)+2)*sqrt(sn(n))) ...
       * hypergeom2F1(2, 1/2, lambda(n)+2, s/sn(n)) * theta(sn(n) - s);
end

if (theta(s - sn(n)) > 0 && abs(sn(n)/s) < (1 + eps(s))) % Check sanity
out = out + (4/3)*(gamma(lambda(n)+3/2)) / gamma(lambda(n)) ...
       * (sn(n)^(3/2)/s^2)...
       * hypergeom2F1(1-lambda(n), 2, 5/2, sn(n)/s) * theta(s - sn(n));
end

end