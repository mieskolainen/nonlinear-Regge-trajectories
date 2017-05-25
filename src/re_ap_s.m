% Baryonic regge trajectory, real part of the slope, Re alpha'(s)
% ------------------------------------------------------------------------
%
%     R. Fiore, L. L. Jenkovszky, F. Paccanoni, A. Prokudin,
%     "Baryonic Regge trajectories with analyticity constraints",
%     Phys.Rev. D70 (2004) 054003
%
% mikael.mieskolainen@cern.ch

function out = re_ap_s(s)

global c;

% Total number of resonances
N = length(c);

out = 0;
for n = 1:N
   out = out + c(n)*B(s,n);
end
out = out / pi;

% Take real part due to floating points
out = real(out);

end

% Amplitudes
function out = B(s,n)

global sn;
global lambda;
global delta;

out = 0;

if (theta(sn(n) - s) > 0 && abs(s/sn(n)) < (1+eps(s)) ) % Check sanity
out = out + (gamma(1-delta)*gamma(lambda(n)+1)) / (gamma(lambda(n)-delta+2)*sn(n)^(1-delta)) ...
      * hypergeom2F1(2, 1-delta, lambda(n)-delta+2, s/sn(n)) * theta(sn(n) - s);
end

if (theta(s - sn(n)) > 0 && abs(sn(n)/s) < (1+eps(s)) ) % Check sanity
out = out + ( pi*s^(delta-1)*((s-sn(n))/s)^lambda(n) * cot(pi*(1-delta)) * (delta + lambda(n)*sn(n)/(s-sn(n))) ...
      - (gamma(-delta)*gamma(lambda(n)+1)*sn(n)^(delta+1)) / ((1+delta)*gamma(lambda(n)-delta)*s^2) ...
      * hypergeom2F1(1+delta-lambda(n), 2, delta+2, sn(n)/s) ) * theta(s - sn(n));
end

end