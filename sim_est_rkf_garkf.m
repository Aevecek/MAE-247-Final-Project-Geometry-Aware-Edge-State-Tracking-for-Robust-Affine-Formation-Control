% sim_est_rkf_garkf.m
clear; clc;

% simulation parameters
dt = 0.01;
t_max = 30;
K_max = floor(t_max / dt);
MC_RUN = 10; % Monte Carlo iterations

% Create directory for results
if ~exist('results', 'dir')
    mkdir('results');
end

% Lambdas across different reliabilities
sense_mat_1  = binornd(1, 1.0, [10, 10, K_max]);
sense_mat_04 = binornd(1, 0.4, [10, 10, K_max]);
sense_mat_02 = binornd(1, 0.2, [10, 10, K_max]);
sense_mat_00 = binornd(1, 1.0, [10, 10, K_max]);
sense_mat_00(2, 1, :) = 0; % Node 1 to Node 2 missing

sim_params.dt = dt;
sim_params.t_max = t_max;
sim_params.static = false;
sim_params.shape = 'hexa';
sim_params.sigma_v = 0.1;

%% RKF
% standard RKF
sim_params.RAL = false;
sim_params.conRAL = false;
sim_params.GARKF = false;
sim_params.RKFIO = true;

lambdas = {1, 0.4, 0.2, 0};
sense_mats = {sense_mat_1, sense_mat_04, sense_mat_02, sense_mat_00};
names = {'1', '04', '02', '00'};

for i = 1:length(lambdas)
    sim_params.sense_mat = sense_mats{i};
    
    rkf_est_track_mc = zeros(2, K_max, MC_RUN);
    true_edge_track_mc = zeros(2, K_max, MC_RUN);
    rkf_cov_track_mc = zeros(6, 6, K_max, MC_RUN);
    
    for j = 1:MC_RUN
        fw = Framework_EST_OUT_OF_LOOP(sim_params);
        fw.run();
        rkf_est_track_mc(:, :, j) = fw.rkf_est_track;
        true_edge_track_mc(:, :, j) = fw.true_edge_track;
        rkf_cov_track_mc(:, :, :, j) = fw.cov_track;
    end
    
    est_error = mean(vecnorm(rkf_est_track_mc - true_edge_track_mc, 2, 1), 3); % scalar mean tracking error
    est_cov = mean(arrayfun(@(k) trace(rkf_cov_track_mc(:,:,k)), 1:K_max), 3); % trace of covariance matrix
    
    save(sprintf('results/est_error_rkf_lmd_%s.mat', names{i}), 'est_error');
    save(sprintf('results/est_cov_rkf_lmd_%s.mat', names{i}), 'est_cov');
end

%% GA-RKF
% for GA-RKF
sim_params.RAL = true;
sim_params.conRAL = true;
sim_params.GARKF = true;
sim_params.RKFIO = false;

for i = 1:length(lambdas)
    sim_params.sense_mat = sense_mats{i};
    
    garkf_est_track_mc = zeros(2, K_max, MC_RUN);
    true_edge_track_mc = zeros(2, K_max, MC_RUN);
    garkf_cov_track_mc = zeros(6, 6, K_max, MC_RUN);
    
    for j = 1:MC_RUN
        fw = Framework_EST_OUT_OF_LOOP(sim_params);
        fw.run();
        garkf_est_track_mc(:, :, j) = fw.con_garkf_est_track;
        true_edge_track_mc(:, :, j) = fw.true_edge_track;
        garkf_cov_track_mc(:, :, :, j) = fw.cov_track;
    end
    
    est_error = mean(vecnorm(garkf_est_track_mc - true_edge_track_mc, 2, 1), 3);
    est_cov = mean(arrayfun(@(k) trace(garkf_cov_track_mc(:,:,k)), 1:K_max), 3);
    
    save(sprintf('results/est_error_garkf_lmd_%s.mat', names{i}), 'est_error'); % scalar mean tracking error
    save(sprintf('results/est_cov_garkf_lmd_%s.mat', names{i}), 'est_cov'); % trace of covariance matrix
end

disp('Simulation done');
