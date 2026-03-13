% sim_grid_blackout.m
clear; clc;

dt = 0.01;
t_max = 60;
K_max = floor(t_max / dt);

% Create directory for results
if ~exist('results', 'dir')
    mkdir('results');
end

% initialise parameters
sim_params.dt = dt;
sim_params.t_max = t_max;
sim_params.static = false;
sim_params.sigma_v = 0.1;
sim_params.shape = 'grid';
sim_params.sense_mat = binornd(1, 1.0, [16, 16, K_max]); 

% GA-RKF
sim_params.RAL = true;
sim_params.conRAL = true; 
sim_params.GARKF = true;
sim_params.RKFIO = false;

sim_params.blind_start = floor(20 / dt); % start of blind window
sim_params.blind_end = floor(35 / dt); % end of blind window
sim_params.blind_node = 16; 

baseline_params = sim_params; % baseline
baseline_params.blind_start = K_max + 1;
baseline_params.blind_end = K_max + 2;
fw_base = Framework_Leader_Blackout(baseline_params);
fw_base.run();
TE_base = fw_base.tracking_error();

fw_blind = Framework_Leader_Blackout(sim_params); fw_blind.run(); % sim blind leader
TE_blind = fw_blind.tracking_error();
traj_blind = fw_blind.config_track;
P_nom = fw_blind.P; % save nominal coordinates

save('results/grid_blind_data.mat', 'TE_base', 'TE_blind', 'traj_blind', 'P_nom', 'sim_params');
disp('Simulation done');
