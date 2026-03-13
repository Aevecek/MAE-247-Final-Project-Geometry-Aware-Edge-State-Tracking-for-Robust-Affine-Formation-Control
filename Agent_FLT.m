classdef Agent_FLT < Agent
    properties
        Phi_now
    end
    
    methods
        function obj = Agent_FLT(ID, sense_mat, P, Z, dt, B_nom, B_und, L, D, is_leader, Z_target, K_max, static)
            obj@Agent(ID, sense_mat, P, Z, dt, B_nom, B_und, L, D, is_leader, Z_target, K_max, static);
        end
        
        function store_Phi(obj, Phi)
            obj.Phi_now = Phi;
        end
        
        function valid_n = valid_neighbors(obj, k)
            valid_neighbor = find(obj.sense_mat(:, obj.ID, k));
            valid_n = intersect(obj.neighbors, valid_neighbor);
        end
    end
end
