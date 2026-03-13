% plot_target_traj.m
clear; clc; close all;

dt = 0.01;
t_max = 60;
K_max = floor(t_max / dt);
sim_params.dt = dt;

sim_params.t_max = t_max;
sim_params.static = false;
sim_params.shape = 'hexa';
sim_params.sigma_v = 0.1;
sim_params.sense_mat = binornd(1, 1, [10, 10, K_max]);

fw = Framework(sim_params);
ref_traj = fw.Z_target;
B = fw.B_nom;
[n, m] = size(B);
edges = zeros(2, m);
for i = 1:m
    idx = find(B(:, i)); edges(:, i) = idx;
end

color_follower = [146, 208, 80] / 256.0; 
color_leader = [255, 184, 28] / 256.0;
figure('Position', [100, 100, 900, 300]);
hold on;

for i = 1:n
    x = squeeze(ref_traj(1, i, :)); 
    y = squeeze(ref_traj(2, i, :));
    if ismember(i, [1, 2, 3])
        plot(x, y, '-', 'Color', color_leader, 'LineWidth', 1);
    else
        plot(x, y, '-', 'Color', color_follower, 'LineWidth', 1);
    end
end

stamps = [2, 1600, 2800, 4000, 6000]; time_labels = [0, 16, 28, 40, 60]; 
for idx_stamp = 1:length(stamps)
    t = stamps(idx_stamp); 
    nodes_xy = squeeze(ref_traj(:, :, t))'; 
    for j = 1:m
        idx = edges(:, j); 
        plot(nodes_xy(idx, 1), nodes_xy(idx, 2), 'k-', 'LineWidth', 1);
    end
    for i = 1:n
        if ismember(i, [4, 5, 7])
            plot(nodes_xy(i, 1), nodes_xy(i, 2), '.', 'Color', color_leader, 'MarkerSize', 20);
        else
            plot(nodes_xy(i, 1), nodes_xy(i, 2), '.', 'Color', color_follower, 'MarkerSize', 20);
        end
    end
    text(mean(nodes_xy(:, 1)), min(nodes_xy(:, 2)) - 0.8, sprintf('$t = %ds$', time_labels(idx_stamp)), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Interpreter', 'latex');
end

xlabel('$x[m]$', 'Interpreter', 'latex', 'FontSize', 14); 
ylabel('$y[m]$', 'Interpreter', 'latex', 'FontSize', 14);
xlim([-5, 48]);
ylim([-5, 5]); 
box on; 
tightfig;

% Create directory for results
if ~exist('results', 'dir')
    mkdir('results');
end
saveas(gcf, 'results/target_traj.png');