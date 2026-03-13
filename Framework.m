classdef Framework < handle
    properties
        dt, K_max, static, shape, B_und, B_nom, N, M_und, M_nom, L
        leader_ID, D, P, Z_init, Z_target, Theta_target, config_track
        U_last, R, obs_noise, sense_mat, agents
    end
    
    methods
        function obj = Framework(sim_params)
            obj.dt = sim_params.dt;
            obj.K_max = floor(sim_params.t_max / sim_params.dt);
            obj.static = sim_params.static;
            obj.shape = sim_params.shape;
            
            [obj.B_und, obj.B_nom, obj.N, obj.D, obj.P] = graph_config(obj.shape);
            obj.M_und = size(obj.B_und, 2);
            obj.M_nom = size(obj.B_nom, 2);
            
            if strcmp(obj.shape, 'hexa')
                obj.leader_ID = [0, 0, 0, 1, 1, 0, 1, 0, 0, 0];
                data = load('data/stress_hexa.mat');
                obj.L = data.L;
            elseif strcmp(obj.shape, 'grid')
                obj.leader_ID = zeros(1, 16);
                obj.leader_ID([1, 4, 13, 16]) = 1; 
                P_aug = [ones(16, 1), obj.P'];
                [U, ~, ~] = svd(P_aug);
                U_null = U(:, 4:end); 
                obj.L = U_null * U_null';
            end
            
            [obj.Z_init, obj.Z_target, obj.Theta_target] = obj.config_setup();
            obj.config_track = zeros(obj.D, obj.N, obj.K_max + 1);
            obj.config_track(:, :, 1) = obj.Z_init;
            obj.U_last = zeros(obj.D, obj.N);
            obj.R = (sim_params.sigma_v^2) * [1, 0.3; 0.3, 1];
            
            noise_flat = mvnrnd([0, 0], obj.R, obj.K_max * obj.N * obj.N);
            obj.obs_noise = reshape(noise_flat, [obj.K_max, obj.N, obj.N, obj.D]);
            obj.sense_mat = sim_params.sense_mat;
            
            obj.agents = cell(1, obj.N);
            for i = 1:obj.N
                obj.agents{i} = Agent(i, obj.sense_mat, obj.P, obj.Z_init, obj.dt, ...
                    obj.B_nom, obj.B_und, obj.L, obj.D, obj.leader_ID(i), ...
                    obj.Z_target, obj.K_max, obj.static);
            end
        end
        
        function [Z_init, Z_target, Theta_target] = config_setup(obj)
            Z_init = mvnrnd(zeros(1, obj.D), 2 * eye(obj.D), obj.N)'; 
            Z_target = zeros(obj.D, obj.N, obj.K_max);
            
            if obj.static
                Theta_target = repmat(eye(2), [1, 1, obj.K_max]);
                target_config = eye(2) * obj.P;
                Z_target = repmat(target_config, [1, 1, obj.K_max]);
            else
                dataA = load('data/true_param.mat'); dataT = load('data/trans.mat');
                step_idx = int32(obj.dt / 0.001);
                A_full = dataA.A; t_full = dataT.t;
                
                Theta_target = A_full(:, :, 1:step_idx:end);
                Theta_target = Theta_target(:, :, 1:obj.K_max);
                trans_target = t_full(:, 1:step_idx:end);
                trans_target = trans_target(:, 1:obj.K_max);
                
                for k = 1:obj.K_max
                    Z_target(:, :, k) = Theta_target(:, :, k) * obj.P + repmat(trans_target(:, k), 1, obj.N);
                end
            end
        end
        
        function error = tracking_error(obj)
            err_matrix = obj.config_track(:, :, 1:end-1) - obj.Z_target;
            error = squeeze(sum(sum(err_matrix.^2, 1), 2)) / obj.N;
        end
        
        function run(obj)
            for k = 1:obj.K_max
                for i = 1:obj.N
                    Zij = obj.config_track(:, i, k) - obj.config_track(:, :, k);
                    Z = obj.config_track(:, :, k);
                    noise_ij = squeeze(obj.obs_noise(k, i, :, :))';
                    Yij = Zij + noise_ij;
                    sel_mat = repmat(obj.sense_mat(:, i, k)', obj.D, 1);
                    Yij = Yij .* sel_mat;
                    [obj.config_track(:, i, k + 1), obj.U_last(:, i)] = ...
                        obj.agents{i}.step(Z, obj.U_last, Yij, k);
                end
            end
        end
    end
end