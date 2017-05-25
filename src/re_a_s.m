% Baryonic trajectory, real part, Re alpha(s)
% ------------------------------------------------------------------------
%
%     R. Fiore, L. L. Jenkovszky, F. Paccanoni, A. Prokudin,
%     "Baryonic Regge trajectories with analyticity constraints",
%     Phys.Rev. D70 (2004) 054003
%
%
% mikael.mieskolainen@cern.ch

function out = re_a_s(s)

global c;
global alpha0;

% Total number of resonances
N = length(c);

out = 0;
for n = 1:N
   out = out + c(n)*A(s,n);
end
out = out * s / pi + alpha0;

% Take real part due to floating points
out = real(out);

end

% Amplitudes
function out = A(s,n)

global sn;
global lambda;
global delta;

out = 0;

cut = 9.999; % Parameter

if (theta(sn(n) - s) > 0 && abs(s/sn(n)) < (1+eps(s))) % Check sanity
    
    %fprintf('1:: %0.3f \n', lambda(n)-delta+2);
out = out + (gamma(1-delta)*gamma(lambda(n)+1)) / (gamma(lambda(n)-delta+2)*sn(n)^(1-delta)) ...
      * hypergeom2F1(1, 1-delta, min(lambda(n)-delta+2, cut), s/sn(n)) * theta(sn(n) - s);
end

if (theta(s - sn(n)) > 0 && abs(sn(n)/s) < (1+eps(s))) % Check sanity
    
    %fprintf('2:: %0.3f \n', delta-lambda(n));
    
out = out + ( pi*s^(delta-1)*((s-sn(n))/s)^lambda(n) * cot(pi*(1-delta)) ...
      - (gamma(-delta)*gamma(lambda(n)+1)*sn(n)^delta) / (s*gamma(lambda(n)-delta+1)) ...
      * hypergeom2F1(max(delta-lambda(n), -cut), 1, delta+1, sn(n)/s) * theta(s - sn(n)) );
end

end
