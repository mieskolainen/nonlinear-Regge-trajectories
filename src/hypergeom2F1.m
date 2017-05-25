% Gauss hypergeometric 2F1(a,b,c,z, precision) by power series expansion
% ------------------------------------------------------------------------
%
% https://en.wikipedia.org/wiki/Hypergeometric_function
% 
% Matlab has also a function hypergeom() for the generalized case.
%
%
% Cross-Checked against: http://keisan.casio.com/exec/system/1349143084
% Set the precision argument, if extra accuracy needed (default 1e-11)
%
% mikael.mieskolainen@cern.ch

function f = hypergeom2F1(a,b,c,z, precision)

if (abs(z) > 1)
    fprintf('hfpergeometric2F1: for z = %0.3f must hold |z| < 1 \n', z);
    return;
end

if (c <= 0 && c == floor(c))
    fprintf('hfpergeometric2F1: c = %0.3f must not be from the singular set {...,-4,-3,-2,-1,0} \n', c)
    return;
end

% A.) First special cases (see Wikipedia):
fs = 5; % Floating point safe scale, semi-arbitrarf

% Gauss's theorem (z = 1)
if (abs(z-1) < eps(z)*fs)
    f = (gamma(c)*gamma(c-a-b)) / (gamma(c-a)*gamma(c-b));
    return;
end

% Kummer's theorem (z = -1)
if (abs(z-(-1)) < eps(z)*fs && abs(c-(1+a-b)) < eps(c)*fs)
    f = (gamma(1 + a - b)*gamma(1 + 0.5*a)) / (gamma(1+a)*gamma(1+0.5*a-b));
    return;
end

% B.) Generic case by truncated series expansion
if (nargin == 4)
    precision = 1e-11; % Default precision
end
delta = 1;
f = 1;
m = 0;

max_order = 1000; % Safety stop

while (precision < abs(delta))
   delta = delta*z*(a + m)*(b + m) / (c + m) / (m + 1);
   f = f + delta;
   m = m + 1;
   
   if (m > max_order)
      break; 
   end
end

end