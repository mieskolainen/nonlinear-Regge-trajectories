% Plot the model decaywidth and DATA
%
% mikael.mieskolainen@cern.ch

function plotwidth(M2, N_M, N_M_err, N_W, N_W_err, names, type, traj)

% Mass^2-range
width = zeros(size(M2));

% Baryonic
if (strcmp(type, 'fermionic'))
    for i = 1:length(M2)
        alpha_prime_real  = re_ap_s(M2(i));
        alpha_imag        = im_a_s(M2(i));
        width(i) = alpha_imag / (sqrt(M2(i))*alpha_prime_real);
    end
% Mesonic
elseif (strcmp(type, 'bosonic'))
    for i = 1:length(M2)
        alpha_prime_real  = re_apM_s(M2(i));
        alpha_imag        = im_aM_s(M2(i));
        width(i) = alpha_imag / (sqrt(M2(i))*alpha_prime_real);
    end
end

% Plot model
plot(M2, width); hold on;

% Plot Data
errorbar(N_M.^2, N_W, N_W_err, 'r.');
errorbar(N_M.^2, N_W', (N_M_err ./ N_M) .* (N_M.^2), 'r.', 'horizontal');
for i = 1:length(N_M)
   text(N_M(i).^2, N_W(i) - 0.01, sprintf('  %s', names{i}), 'fontsize', 8, 'rotation', -20);
end
xlabel('$M^2$ (GeV$^2$)', 'interpreter', 'latex');
ylabel('$\Gamma = $ Im $\alpha(M^2) / (M$Re $\alpha''(M^2))$ (GeV)','interpreter','latex');

title(sprintf('$%s$-trajectory: width', traj),'interpreter','latex');
axis([0 M2(end) 0 0.7]);
set(gca,'fontsize', 12);

end