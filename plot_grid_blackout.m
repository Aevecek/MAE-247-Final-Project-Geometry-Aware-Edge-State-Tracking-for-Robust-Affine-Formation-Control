% plot_grid_blackout.m
clear; clc; close all;

load('results/grid_blind_data.mat');
dt = sim_params.dt;
t_max = sim_params.t_max;
K_max = floor(t_max / dt);
t_array = linspace(0, t_max, K_max)';

blind_time_start = sim_params.blind_start * dt;
blind_time_end = sim_params.blind_end * dt;

color_follower = [146, 208, 80] / 256.0;
color_leader = [255, 184, 28] / 256.0;
color_blind = [255, 50, 50] / 256.0;

figure('Position', [100, 100, 1000, 400]);
y_max = max(max(TE_base), max(TE_blind)) * 5; 
y_min = max(1e-6, min(min(TE_base), min(TE_blind)) * 0.1); 

subplot(1, 2, 1); hold on;
fill([blind_time_start, blind_time_end, blind_time_end, blind_time_start], ...
     [y_min, y_min, y_max, y_max], [0.9 0.9 0.9], 'EdgeColor', 'none', 'DisplayName', 'Leader 16 GPS Failure');
plot(t_array, TE_base, 'Color', [0.6 0.8 0.2], 'LineWidth', 2, 'DisplayName', 'Baseline (Perfect)');
plot(t_array, TE_blind, 'Color', [0.86 0.08 0.24], 'LineWidth', 2, 'DisplayName', 'Top-Right Leader Blindfolded');
set(gca, 'YScale', 'log');
ylim([y_min, y_max]);
xlim([0, 60]); 
title('4x4 Grid Tracking Error (\delta_k)');
xlabel('Time (s)'); ylabel('Error');
legend('Location', 'northeast');
grid on;
box on;

subplot(1, 2, 2); hold on;
for i = 1:16
    x = squeeze(traj_blind(1, i, :));
    y = squeeze(traj_blind(2, i, :));
    if i == sim_params.blind_node
        plot(x, y, '-', 'Color', color_blind, 'LineWidth', 1.5);
    elseif ismember(i, [1, 4, 13])
        plot(x, y, '-', 'Color', color_leader, 'LineWidth', 1);
    else
        plot(x, y, '-', 'Color', color_follower, 'LineWidth', 1);
    end
end

stamps = [2, 1600, 2800, 4000, 6000];
time_labels = [0, 16, 28, 40, 60];
for idx_stamp = 1:length(stamps)
    t = stamps(idx_stamp);
    nodes_xy = squeeze(traj_blind(:, :, t))'; 
    for i = 1:16
        for j = i+1:16
            if norm(P_nom(:,i) - P_nom(:,j)) == 1
                plot([nodes_xy(i,1), nodes_xy(j,1)], [nodes_xy(i,2), nodes_xy(j,2)], 'k-', 'LineWidth', 1.5);
            end
        end
    end
    for i = 1:16
        if i == sim_params.blind_node
            plot(nodes_xy(i, 1), nodes_xy(i, 2), '.', 'Color', color_blind, 'MarkerSize', 20);
        elseif ismember(i, [1, 4, 13])
            plot(nodes_xy(i, 1), nodes_xy(i, 2), '.', 'Color', color_leader, 'MarkerSize', 15);
        else
            plot(nodes_xy(i, 1), nodes_xy(i, 2), '.', 'Color', color_follower, 'MarkerSize', 15);
        end
    end
    center_x = mean(nodes_xy(:, 1));
    min_y = min(nodes_xy(:, 2));
    if time_labels(idx_stamp) == 28
        text(center_x, min_y - 1.5, sprintf('$t = %ds$ (GPS Lost)', time_labels(idx_stamp)), 'HorizontalAlignment', 'center', 'FontSize', 10, 'Interpreter', 'latex', 'Color', 'r');
    else
        text(center_x, min_y - 1.5, sprintf('$t = %ds$', time_labels(idx_stamp)), 'HorizontalAlignment', 'center', 'FontSize', 10, 'Interpreter', 'latex');
    end
end
title('4x4 Grid Trajectory Progression');
xlabel('$x[m]$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y[m]$', 'Interpreter', 'latex', 'FontSize', 12);
xlim([-5, 48]);
ylim([-8, 8]);
box on;
tightfig;
saveas(gcf, 'results/grid_blind_test.png');