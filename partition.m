function B_bidir = partition(B_bid, i, shape)
    if strcmp(shape, 'hexa')
        Ms = [5, 5, 5, 8, 7, 7, 8, 5, 5, 5]; % Hexagon shape
    elseif strcmp(shape, 'grid')
        Ms = repmat(15, 1, 16); % Grid shape
    else
        error('Unknown shape');
    end

    idx_start = sum(Ms(1:i-1)) + 1;
    idx_end = idx_start + Ms(i) - 1;
    B_bidir = B_bid(:, idx_start:idx_end);
end