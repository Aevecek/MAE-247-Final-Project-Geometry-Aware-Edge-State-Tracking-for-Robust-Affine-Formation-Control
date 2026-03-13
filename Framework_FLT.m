classdef Framework_FLT < Framework
    properties
        Theta_est_track
        Theta_est_now
        Theta_conRAL_track
        use_RAL
        use_conRAL
        use_GARKF
        use_RKFIO
        Phi_track
        psi_track
        ci_track
        bi_track
        sigma_w
        state
        cov
        cov_track
    end
    
    methods
        function obj = Framework_FLT(sim_params)
            obj@Framework(sim_params);
            
            obj.Theta_est_track = zeros(obj.D, obj.D, obj.N, obj.K_max);
            obj.Theta_est_now = zeros(obj.D, obj.D, obj.N);
            obj.Theta_conRAL_track = zeros(obj.D, obj.D, obj.N, obj.K_max + 1);
            
            for i = 1:obj.N
                obj.agents{i} = Agent_FLT(i, obj.sense_mat, obj.P, obj.Z_init, obj.dt, ...
                    obj.B_nom, obj.B_und, obj.L, obj.D, obj.leader_ID(i), ...
                    obj.Z_target, obj.K_max, obj.static);
            end
            
            obj.use_RAL = sim_params.RAL;
            obj.use_conRAL = sim_params.conRAL;
            obj.use_GARKF = sim_params.GARKF;
            obj.use_RKFIO = sim_params.RKFIO;
            
            obj.Phi_track = cell(1, obj.N);
            for i = 1:obj.N
                obj.Phi_track{i} = zeros(obj.D, length(obj.agents{i}.neighbors), obj.K_max);
            end
            
            obj.psi_track = zeros(obj.N, obj.K_max);
            obj.ci_track = zeros(obj.N, obj.K_max);
            obj.bi_track = zeros(obj.N, obj.K_max);
            
            obj.sigma_w = 0.01;
            obj.state = zeros(obj.N, obj.N, 6);
            obj.cov = repmat(2 * eye(6), [1, 1, obj.N, obj.N]); % [6, 6, N, N]
            obj.cov = permute(obj.cov, [3, 4, 1, 2]); % [N, N, 6, 6]
            
            obj.cov_track = zeros(6, 6, obj.K_max);
        end
        
        function est_Theta(obj, k)
            for i = 1:obj.N
                Bi = partition(obj.B_nom, i, obj.shape);
                sel_seq = obj.sense_mat(obj.agents{i}.neighbors, i, k);
                idx_sel = find(sel_seq);
                Bik = Bi(:, idx_sel);
                
                Zij = obj.config_track(:, i, k) - obj.config_track(:, :, k);
                noise_ij = squeeze(obj.obs_noise(k, i, :, :))';
                Yij = Zij + noise_ij;
                
                Yij_neib = Yij(:, obj.agents{i}.neighbors);
                Yij_sel = Yij_neib(:, idx_sel);
                
                Hi = obj.P * Bik;
                
                if rank(Hi * Hi') == obj.D
                    Phi = (Hi * Hi') \ Hi;
                    Theta = Yij_sel * Phi';
                    obj.Theta_est_track(:, :, i, k) = Theta;
                    obj.Theta_est_now(:, :, i) = Theta;
                    obj.agents{i}.store_Phi(Phi);
                else
                    if k > 1
                        obj.Theta_est_now(:, :, i) = obj.Theta_conRAL_track(:, :, i, k - 1);
                    else
                        obj.Theta_est_now(:, :, i) = zeros(obj.D, obj.D);
                    end
                end
            end
        end
        
        function Zij_est = RAL(obj, i, k, Yij)
            sel_seq = obj.sense_mat(:, i, k);
            idx_sel_seq = find(sel_seq);
            
            if obj.use_conRAL
                Theta = obj.Theta_conRAL_track(:, :, i, k);
            else
                Theta = obj.Theta_est_track(:, :, i, k);
            end
            
            Zij_est = Theta * (obj.P(:, i) - obj.P);
            Zij_est(:, idx_sel_seq) = Yij(:, idx_sel_seq); % Data consistency
        end
        
        function CI(obj, k)
            for i = 1:obj.N
                Nik = length(obj.agents{i}.neighbors);
                psii = 0;
                ci = 0;
                bi = 0;
                
                for j = obj.agents{i}.neighbors
                    psii = psii + norm(obj.Theta_est_track(:, :, i, k) - obj.Theta_est_track(:, :, j, k), 'fro')^2;
                end
                
                obj.bi_track(i, k) = (trace(obj.R) / Nik) * bi;
                obj.ci_track(i, k) = (obj.N / Nik) * ci;
                obj.psi_track(i, k) = (1 / Nik) * psii;
            end
        end
        
        function ConRAL(obj, k)
            alp = 0.08;
            if k == 1
                obj.Theta_conRAL_track(:, :, :, 1) = obj.Theta_est_track(:, :, :, 1);
            end
            
            for i = 1:obj.N
                sel_seq = obj.sense_mat(:, i, k);
                idx_sel_seq = find(sel_seq);
                
                neib_sum = zeros(obj.D, obj.D);
                neib_sum1 = zeros(obj.D, obj.D);
                
                for j = obj.agents{i}.neighbors
                    if ismember(j, idx_sel_seq)
                        neib_sum = neib_sum + obj.Theta_conRAL_track(:, :, j, k) - obj.Theta_conRAL_track(:, :, i, k);
                        if any(abs(obj.Theta_est_track(:, :, j, k)) > 1e-5, 'all')
                            neib_sum1 = neib_sum1 + obj.Theta_est_track(:, :, j, k) - obj.Theta_conRAL_track(:, :, i, k);
                        end
                    end
                end
                
                if any(abs(obj.Theta_est_track(:, :, i, k)) > 1e-5, 'all')
                    neib_sum1 = neib_sum1 + obj.Theta_est_track(:, :, i, k) - obj.Theta_conRAL_track(:, :, i, k);
                end
                
                obj.Theta_conRAL_track(:, :, i, k + 1) = obj.Theta_conRAL_track(:, :, i, k) + alp * (neib_sum + neib_sum1);
            end
        end
        
        function [cstate, ccov, meas] = KF(obj, state, cov, yij, R_val, is_obs)
            dt_val = obj.dt;
            Q = (obj.sigma_w^2) * [0.25*dt_val^4, 0.5*dt_val^3, 0.5*dt_val^2, 0, 0, 0; ...
                                   0.5*dt_val^3, dt_val^2, dt_val, 0, 0, 0; ...
                                   0.5*dt_val^2, dt_val, 1, 0, 0, 0; ...
                                   0, 0, 0, 0.25*dt_val^4, 0.5*dt_val^3, 0.5*dt_val^2; ...
                                   0, 0, 0, 0.5*dt_val^3, dt_val^2, dt_val; ...
                                   0, 0, 0, 0.5*dt_val^2, dt_val, 1];
                                   
            G = [1, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0];
            F = [1, dt_val, 0.5*dt_val^2, 0, 0, 0; ...
                 0, 1, dt_val, 0, 0, 0; ...
                 0, 0, 1, 0, 0, 0; ...
                 0, 0, 0, 1, dt_val, 0.5*dt_val^2; ...
                 0, 0, 0, 0, 1, dt_val; ...
                 0, 0, 0, 0, 0, 1];
                 
            pstate = F * state(:);
            pcov = F * squeeze(cov) * F' + Q;
            
            if is_obs
                K = pcov * G' / (R_val + G * pcov * G');
                cstate = pstate + K * (yij(:) - G * pstate);
                ccov = (eye(6) - K * G) * pcov;
            else
                cstate = pstate;
                ccov = pcov;
            end
            
            meas = G * cstate;
        end
        
        function Zij_est = GARKF(obj, i, k, Yij)
            sel_seq = obj.sense_mat(:, i, k);
            idx_sel_seq = find(sel_seq);
            Zij_est = zeros(obj.D, obj.N);
            
            for j = obj.agents{i}.neighbors
                if ismember(j, idx_sel_seq)
                    R_val = obj.R;
                else
                    R_val = obj.R + obj.psi_track(i, k) * eye(obj.D);
                    if obj.use_conRAL
                        R_val = (1 / obj.N) * R_val;
                    end
                end
                
                [cst, ccv, meas] = obj.KF(squeeze(obj.state(i, j, :)), squeeze(obj.cov(i, j, :, :)), Yij(:, j), R_val, true);
                obj.state(i, j, :) = cst;
                obj.cov(i, j, :, :) = ccv;
                Zij_est(:, j) = meas;
                
                if i == 1 && j == 2
                    obj.cov_track(:, :, k) = ccv;
                end
            end
        end
        
        function Zij_est = RKFIO(obj, i, k, Yij)
            sel_seq = obj.sense_mat(:, i, k);
            idx_sel_seq = find(sel_seq);
            Zij_est = zeros(obj.D, obj.N);
            
            for j = obj.agents{i}.neighbors
                is_obs = ismember(j, idx_sel_seq);
                [cst, ccv, meas] = obj.KF(squeeze(obj.state(i, j, :)), squeeze(obj.cov(i, j, :, :)), Yij(:, j), obj.R, is_obs);
                obj.state(i, j, :) = cst;
                obj.cov(i, j, :, :) = ccv;
                Zij_est(:, j) = meas;
                
                if i == 1 && j == 2
                    obj.cov_track(:, :, k) = ccv;
                end
            end
        end
        
        function run(obj)
            for k = 1:obj.K_max
                obj.est_Theta(k);
                obj.CI(k);
                if obj.use_conRAL
                    obj.ConRAL(k);
                end
                
                for i = 1:obj.N
                    Zij = obj.config_track(:, i, k) - obj.config_track(:, :, k);
                    Z = obj.config_track(:, :, k);
                    
                    noise_ij = squeeze(obj.obs_noise(k, i, :, :))';
                    Yij = Zij + noise_ij;
                    
                    sel_mat = repmat(obj.sense_mat(:, i, k)', obj.D, 1);
                    Yij = Yij .* sel_mat;
                    
                    if obj.use_RAL
                        Zij_est = obj.RAL(i, k, Yij);
                        if obj.use_GARKF
                            Zij_est = obj.GARKF(i, k, Zij_est);
                        end
                    elseif obj.use_RKFIO
                        Zij_est = obj.RKFIO(i, k, Yij);
                    else
                        Zij_est = Yij;
                    end
                    
                    [obj.config_track(:, i, k + 1), obj.U_last(:, i)] = ...
                        obj.agents{i}.step(Z, obj.U_last, Zij_est, k);
                end
            end
        end
    end
end
