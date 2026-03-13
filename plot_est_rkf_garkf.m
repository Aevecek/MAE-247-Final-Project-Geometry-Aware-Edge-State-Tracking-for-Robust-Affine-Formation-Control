% plot_est_rkf_garkf.m
clear; clc; close all;

lambdas = {'1', '04', '02', '00'};
lmd_vals = [1.0, 0.4, 0.2, 0.0];
err_rkf = zeros(1, 4); err_garkf = zeros(1, 4);

for i = 1:length(lambdas)
    data = load(['results/est_error_rkf_lmd_', lambdas{i}, '.mat']);
    err_rkf(i) = mean(data.est_error);
    data = load(['results/est_error_garkf_lmd_', lambdas{i}, '.mat']);
    err_garkf(i) = mean(data.est_error);
end

dt = 0.01;
t_max = 30;
K_max = floor(t_max / dt);
t = linspace(0, t_max, K_max);
cov_lambdas = {'00', '02', '1'};
cov_titles = {'\lambda = 0.0', '\lambda = 0.2', '\lambda = 1.0'};
cov_ylims = {[1e-3, 1e6], [5e-3, 1e1], [1e-2, 1e1]};

figure('Position', [100, 400, 900, 300]); 
for j = 1:length(cov_lambdas)
    lmd_str = cov_lambdas{j};
    load(['results/est_cov_rkf_lmd_', lmd_str, '.mat'], 'est_cov');
    cov_rkf = est_cov;
    load(['results/est_cov_garkf_lmd_', lmd_str, '.mat'], 'est_cov');
    cov_garkf = est_cov;
    
    subplot(1, 3, j);
    plot(t, cov_rkf, '-', 'Color', [0.6 0.8 0.2], 'LineWidth', 1.5, 'DisplayName', 'RKF');
    hold on;
    plot(t, cov_garkf, ':', 'Color', [0.86 0.08 0.24], 'LineWidth', 1.5, 'DisplayName', 'GA-RKF');
    set(gca, 'YScale', 'log'); ylim(cov_ylims{j}); xlim([0, 30]);
    title(cov_titles{j}); xlabel('Time (s)');
    if j == 1 
        ylabel('tr(\Lambda_{ij,k})'); 
    end
    legend('Location', 'best');
    grid on;
end
tightfig;
saveas(gcf, 'results/est_cov_subplots.png');