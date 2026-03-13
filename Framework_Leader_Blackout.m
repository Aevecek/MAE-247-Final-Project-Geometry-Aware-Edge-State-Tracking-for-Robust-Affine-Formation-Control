classdef Framework_Leader_Blackout < Framework_FLT
    properties
        blind_start, blind_end, blind_node
    end
    methods
        function obj = Framework_Leader_Blackout(sim_params)
            obj@Framework_FLT(sim_params);
            obj.blind_start = sim_params.blind_start;
            obj.blind_end = sim_params.blind_end;
            obj.blind_node = sim_params.blind_node; 
        end
        function run(obj)
            for k = 1:obj.K_max
                if k >= obj.blind_start && k <= obj.blind_end
                    obj.agents{obj.blind_node}.is_leader = 0; 
                else
                    obj.agents{obj.blind_node}.is_leader = 1; 
                end
                
                obj.est_Theta(k); obj.CI(k);
                if obj.use_conRAL, obj.ConRAL(k); end
                
                for i = 1:obj.N
                    Zij = obj.config_track(:, i, k) - obj.config_track(:, :, k);
                    Z = obj.config_track(:, :, k);
                    noise_ij = squeeze(obj.obs_noise(k, i, :, :))';
                    Yij = Zij + noise_ij;
                    sel_mat = repmat(obj.sense_mat(:, i, k)', obj.D, 1);
                    Yij = Yij .* sel_mat;
                    
                    if obj.use_RAL
                        Zij_est = obj.RAL(i, k, Yij);
                        if obj.use_GARKF, Zij_est = obj.GARKF(i, k, Zij_est); end
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