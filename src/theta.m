% Heaviside theta function, R -> {0, 0.5, 1}
% ------------------------------------------------------------------------
% Careful comparison of floating point numbers is done.
%
% input       x = input value
%        origin = 0, 0.5 or 1 (default 0.5)
%
% mikael.mieskolainen@cern.ch

function y = theta(x, origin)

if (nargin == 1)
   origin = 0.5; % Default 
end

y = x;

for i = 1:length(x)
    
    if (x(i) > eps(x(i))*4)
        y(i) = 1;
    elseif (x(i) < -eps(x(i))*4)
        y(i) = 0;
    else % x = 0
        y(i) = origin;
    end
end

end