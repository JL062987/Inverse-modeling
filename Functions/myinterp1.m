function [yq] = myinterp1(x, y, xq)
    % myinterp1 1-D periodic interpolation along the 3rd dimension of a 3-D array.
    %
    % Inputs:
    %   x  - Vector of sample points. Must be a subset of [0, 1). Length must match size(y, 3).
    %   y  - 3-D array of sample values (size: n x m x P). Assumes periodic boundary conditions along the 3rd dimension.
    %   xq - Vector of query points. Must be a subset of [0, 1).
    %
    % Outputs:
    %   yq - 3-D array of interpolated values (Size: n x m x length(xq)).

    % Ensure x and xq are column vectors
    x = x(:);
    xq = xq(:);

    % Enforce periodicity by padding the sample points to the left (-1) and right (+1)
    x_padded = [x(end) - 1; x; x(1) + 1];    
    
    % Bring the target 3rd dimension to the front (New size: P x n x m)
    y_perm = permute(y, [3, 1, 2]);

    % Duplicate and tile the y-data to match the padded x-intervals
    y_padded = [y_perm(end,:,:); y_perm; y_perm(1,:,:)];
    
    % Perform highly optimized vectorized 1-D interpolation
    yq_perm = interp1(x_padded, y_padded, xq, 'linear');

    % Restore the original matrix structure (Returned size: n x m x length(xq))
    yq = ipermute(yq_perm, [3, 1, 2]);

end
