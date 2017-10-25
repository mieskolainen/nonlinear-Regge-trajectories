% Plot the Regge trajectory in the complex plane
%
% mikael.mieskolainen@cern.ch, 2017

function plotcomplexplane(s_range, type, traj)

alpha_real = zeros(size(s_range));
alpha_imag = zeros(size(s_range));

% Baryonic
if (strcmp(type,'fermionic'))
    for i = 1:length(s_range)
        alpha_real(i) = re_a_s(s_range(i));
        alpha_imag(i) = im_a_s(s_range(i));
    end

% Mesonic
elseif (strcmp(type,'bosonic'))
   for i = 1:length(s_range)
        alpha_real(i) = re_aM_s(s_range(i));
        alpha_imag(i) = im_aM_s(s_range(i));
   end
end

% Plot model
plot(alpha_real, alpha_imag, 'r-', 'linewidth', 1.25); hold on;

xlabel('Re $\alpha(s)$', 'interpreter', 'latex');
ylabel('Im $\alpha(s)$', 'interpreter','latex');
title(sprintf('$%s$-trajectory', traj), 'interpreter', 'latex');
set(gca,'fontsize', 12);
axis square;
axis tight;
%grid on;

% Draw x-y-axis
x = linspace(s_range(1), s_range(end), 20);
y = linspace(s_range(1), s_range(end)/3, 20); % imaginary part usually small -> / 3
plot(x, zeros(size(x)), 'k-'); hold on;
plot(zeros(size(y)), y, 'k-');

% Draw s values
step = round(length(s_range) / 10);
for i = 1:step:length(s_range)
    plot(alpha_real(i), alpha_imag(i), 'k.'); 
    h = text(1.02*alpha_real(i), 1.1*alpha_imag(i), sprintf('$%0.1f$', s_range(i)), ...
       'fontsize', 9, 'interpreter','latex'); 
    set(h,'Rotation',-45);
end

l = legend('$\alpha(s) = J$', 'location','NE');
legend('boxoff'); set(l,'interpreter','latex');

end