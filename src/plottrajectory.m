% Plot the Regge trajectory and DATA
%
% mikael.mieskolainen@cern.ch

function plottrajectory(s_range, N_M, N_M_err, N_J, alpha0, names, type, traj)

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
plot(s_range, alpha_real, 'k-'); hold on;
plot(s_range, alpha_imag, 'b--'); hold on;

% Plot data
errorbar(N_M.^2, N_J', (N_M_err ./ N_M) .* (N_M.^2), 'r.', 'horizontal');
for i = 1:length(N_M)
   text(N_M(i).^2, N_J(i) + 0.2, sprintf('  %s', names{i}), 'fontsize', 8);
end

xlabel('$s$ (GeV$^2$)', 'interpreter', 'latex');
ylabel('$\alpha(s) = J$','interpreter','latex');
title(sprintf('$%s$-trajectory', traj), 'interpreter', 'latex');
set(gca,'fontsize', 12);
axis square;
axis tight;
%grid on;

% Draw x-y-axis
x = linspace(s_range(1), s_range(end), 20);
y = linspace(s_range(1), s_range(end), 20);
plot(x, zeros(size(x)), 'k-'); hold on;
plot(zeros(size(y)), y, 'k-');

l = legend(sprintf('Re $\\alpha(s), \\alpha_0 = %0.3f$', alpha0), 'Im $\alpha(s)$', 'location','NW');
legend('boxoff'); set(l,'interpreter','latex');

end