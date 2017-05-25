% Baryonic trajectory, imaginary part, Im alpha(s)
% ------------------------------------------------------------------------
%
%     R. Fiore, L. L. Jenkovszky, F. Paccanoni, A. Prokudin,
%     "Baryonic Regge trajectories with analyticity constraints",
%     Phys.Rev. D70 (2004) 054003
%
% mikael.mieskolainen@cern.ch

function out = im_a_s(s)

global sn;
global c;
global lambda;
global delta;

% Total number of resonances
N = length(sn);

out = 0;
for n = 1:N
    if (theta(s - sn(n)) > 0) % Check sanity
        out = out + c(n) * ((s - sn(n))/s)^lambda(n) * theta(s - sn(n));
    end
end

out = s^delta * out;

% Take real part due to floating points
out = real(out);

end