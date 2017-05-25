% Mesonic trajectory, real part, Re alpha(s)
% ------------------------------------------------------------------------
%
%     "Analytic Model of Regge Trajectories"
%     R. Fiore, L. Jenkovszky, V. Magas, F. Paccanoni, A. Papa,
%     Eur.Phys.J. A (2001)
%
% mikael.mieskolainen@cern.ch

function out = re_aM_s(s)

global c;
global alpha0;

% Total number of resonances
N = length(c);

out = 0;
for n = 1:N
   out = out + c(n)*A(s,n);
end
out = out / sqrt(pi) + alpha0;

% Take real part due to floating points
out = real(out);

end

% Amplitudes
function out = A(s,n)

global sn;
global lambda;

out = 0;

if (theta(sn(n) - s) > 0 && abs(s/sn(n)) < (1 + eps(s)) ) % Check sanity
out = out + (s*gamma(lambda(n)+3/2)) / (gamma(lambda(n)+2)*sqrt(sn(n))) ...
       * hypergeom2F1(1, 1/2, lambda(n)+2, s/sn(n)) * theta(sn(n) - s);
end

if (theta(s - sn(n)) > 0 && abs(sn(n)/s) < (1 + eps(s)) ) % Check sanity
out = out + (2*(gamma(lambda(n)+3/2)) / (gamma(lambda(n)+1))*sqrt(sn(n))) ...
       * hypergeom2F1(-lambda(n), 1, 3/2, sn(n)/s) * theta(s - sn(n));
end

end