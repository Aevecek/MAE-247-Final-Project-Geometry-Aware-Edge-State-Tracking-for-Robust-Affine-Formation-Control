function B_bidir = partition(B_bid, i, shape)
% partition extracts the columns from the Bidirectional incidence Matrix
% and extracts the column relevant to agent i
    if strcmp(shape, 'hexa')
        Ms = [5, 5, 5, 8, 7, 7, 8, 5, 5, 5]; % Hexagon shape neighbor counts
    elseif strcmp(shape, 'grid')
        Ms = repmat(15, 1, 16); % Grid shape everyone in communication
    else
        error('Unknown shape'); % error occured
    end

    % slicing the matrix and column indices for agent i
    idx_start = sum(Ms(1:i-1)) + 1;
    idx_end = idx_start + Ms(i) - 1;
    B_bidir = B_bid(:, idx_start:idx_end);
end
