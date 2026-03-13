classdef Agent < handle
    properties
        dt
        sense_mat
        B_nom
        B_und
        N
        M
        L
        ID
        is_leader
        D
        P
        z
        neighbors
        Z_target
        K_max
        static
    end
    
    methods
        function obj = Agent(ID, sense_mat, P, Z, dt, B_nom, B_und, L, D, is_leader, Z_target, K_max, static)
            obj.dt = dt;
            obj.sense_mat = sense_mat;
            obj.B_nom = B_nom;
            obj.B_und = B_und;
            [obj.N, obj.M] = size(obj.B_nom);
            obj.L = L; 
            obj.ID = ID; % 1-based indexing
            obj.is_leader = is_leader;
            obj.D = D;
            obj.P = P;
            obj.z = Z(:, obj.ID); 
            obj.Z_target = Z_target; 
            obj.K_max = K_max;
            obj.static = static;
            
            obj.neighbors = obj.get_neighbors(); 
        end
        
        function neighbor_IDs = get_neighbors(obj)
            % Find edges connected to this agent
            edges = find(obj.B_und(obj.ID, :));
            
            % Find all agents connected to those edges
            [connected_agents, ~] = find(obj.B_und(:, edges));
            
            % Remove duplicates and self
            neighbor_IDs = unique(connected_agents);
            neighbor_IDs(neighbor_IDs == obj.ID) = [];
            
            % Ensure row vector for easy iteration
            neighbor_IDs = neighbor_IDs(:)'; 
        end
        
        function [z_new, u] = step(obj, Y, U, Yij, k)
            if obj.is_leader == 0
                u = zeros(obj.D, 1);
                if obj.static
                    for j = obj.neighbors
                        u = u + 50 * obj.L(obj.ID, j) * Yij(:, j);
                    end
                else
                    gamma = 0;
                    for j = obj.neighbors
                        gamma = gamma + obj.L(obj.ID, j);
                        u = u + obj.L(obj.ID, j) * (Yij(:, j) - U(:, j));
                    end
                    u = -(1 / gamma) * u;
                end
            else
                if obj.static
                    u = -50 * (Y(:, obj.ID) - obj.Z_target(:, obj.ID, k));
                else
                    if k == 1
                        dz_ref = zeros(obj.D, 1);
                    else
                        dz_ref = (obj.Z_target(:, obj.ID, k) - obj.Z_target(:, obj.ID, k - 1)) / obj.dt;
                    end
                    u = -10 * (Y(:, obj.ID) - obj.Z_target(:, obj.ID, k)) + dz_ref;
                end
            end
            
            obj.z = obj.z + obj.dt * u;
            z_new = obj.z;
        end
    end
end