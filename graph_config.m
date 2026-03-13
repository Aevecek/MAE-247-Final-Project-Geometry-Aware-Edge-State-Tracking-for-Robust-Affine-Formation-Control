function [B_und, B_bid, N, D, P] = graph_config(shape)
% graph_config builds the network topology and outputs 
% B_und: Undirected incidence matrix
% B_bid: Bidirectional incidence matrix
% N: Number of nodes
% D: The dimension space
% P: Nominal formation coordinates

    if strcmp(shape, 'hexa') % Load paper's hexagon topology
        B_und = load('data/inc_hexagon.txt');
        B_bid = load('data/inc_hexagon_2way.txt');
        N = 10; D = 2;
        P = [3, 0; 2, 2; 2, -2; 1, 0; ...
            0, 2; 0, -2; -1, 0; -2, 2; ...
            -2, -2; -3, 0]';
    elseif strcmp(shape, 'grid') % Load 4x4 grid topology
        N = 16; D = 2;
        [X, Y] = meshgrid(0:3, 0:3); % generate 4x4 grid
        P = [X(:)'; Y(:)'];
        num_und = 16 * 15 / 2; % build B_und
        B_und = zeros(16, num_und);
        k = 1;
        for i = 1:16
            for j = (i+1):16
                B_und(i, k) = 1; B_und(j, k) = -1; k = k + 1;
            end
        end
        num_nom = 16 * 15; % build B_bid
        B_bid = zeros(16, num_nom);
        k = 1;
        for i = 1:16
            for j = 1:16
                if i ~= j
                    B_bid(i, k) = 1; B_bid(j, k) = -1; k = k + 1;
                end
            end
        end
    else
        error('Unknown Shape'); % error occured
    end
end
