function [M,idx] = Monodromy(X,opt)

    % Compute the monodromy matrix for X.
    % Unless specified, X is assumed to be 1-periodic.

    if nargin < 2
        opt = struct();
    end

    [n,~,N] =  size(X);

    if isfield(opt,'dt')
        dt = opt.dt;
    else
        dt = 1/N;
    end

    
    M = eye(n);
    
    for j = 1:N
        M = expm(X(:,:,j)*dt)*M; % Compute the monodromy matrix
    end
    
    idx = sum( abs( eig(M) ) > 1 ) > 0 ; % Check whether it represents a stable system or not.

end