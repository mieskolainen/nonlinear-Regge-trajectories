% Mesonic trajectory, imaginary part, Im alpha(s)
% ------------------------------------------------------------------------
%
%     "Analytic Model of Regge Trajectories"
%     R. Fiore, L. Jenkovszky, V. Magas, F. Paccanoni, A. Papa,
%     Eur.Phys.J. A (2001)
%
% mikael.mieskolainen@cern.ch

function out = im_aM_s(s)

global sn;
global c;
global lambda;

% Total number of resonances
N = length(c);

out = 0;
for n = 1:N
    if (theta(s - sn(n)) > 0) % Check sanity
    out = out + c(n) * (s - sn(n))^(1/2) ...
        * ((s - sn(n))/s)^lambda(n) * theta(s - sn(n));
    end
end

% Take real part due to floating points
out = real(out);

end