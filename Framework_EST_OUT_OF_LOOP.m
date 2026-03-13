classdef Framework_EST_OUT_OF_LOOP < Framework_FLT
    properties
        rkf_est_track
        con_garkf_est_track
        true_edge_track
    end
    
    methods
        function obj = Framework_EST_OUT_OF_LOOP(sim_params)
            obj@Framework_FLT(sim_params);
            
            obj.rkf_est_track = zeros(obj.D, obj.K_max);
            obj.con_garkf_est_track = zeros(obj.D, obj.K_max);
            obj.true_edge_track = zeros(obj.D, obj.K_max);
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
                   
                    [obj.config_track(:, i, k + 1), obj.U_last(:, i)] = ...
                        obj.agents{i}.step(Z, obj.U_last, Yij, k);
                    
                    if i == 1
                        obj.true_edge_track(:, k) = Zij(:, 2);
                    end
                   
                    sel_mat = repmat(obj.sense_mat(:, i, k)', obj.D, 1);
                    Yij_sel = Yij .* sel_mat;
                    
                    if obj.use_RAL
                        Zij_est = obj.RAL(i, k, Yij_sel);
                        if obj.use_GARKF
                            Zij_est = obj.GARKF(i, k, Zij_est);
                            if i == 1
                                obj.con_garkf_est_track(:, k) = Zij_est(:, 2);
                            end
                        end
                    end
                    
                    if obj.use_RKFIO
                        Zij_est = obj.RKFIO(i, k, Yij_sel);
                        if i == 1
                            obj.rkf_est_track(:, k) = Zij_est(:, 2);
                        end
                    end
                end
            end
        end
    end
end
