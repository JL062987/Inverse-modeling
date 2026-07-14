function [X_fista,h] = CS_nLIM_Group_Lasso_FISTA(A,b,G,lambda,opts)

    [~,N] = size(A);

    % cost function 
    % norm2 = @(v) v'*v; % norm(v)^2
    % calc_F = @(X) 0.5*norm2(A*X-b) + lambda*norm(X,1);
    calc_F = [];
    
    % fista solution 
    opts.pos = false;
    opts.lambda = lambda;

    alpha = 0*lambda; % Ridge coefficient
    
    % graient
    if issparse(A) == 1 && nnz(A)/numel(A) < 0.15
        grad = @(x) A'*(A*x - b) + 2*alpha*x;
    else
        % AtA = full(A'*A);
        % Atb = full(A'*b);
        % grad = @(x) AtA*x - Atb  + 2*alpha*x;
        A = full(A);
        grad = @(x) A'*(A*x - b) + 2*alpha*x;
    end
    
    % initial condition
    Xinit = ones(N,1);

    % Lipschitz constant
    L = svds(A, 1)^2 + 2*alpha;

    if isnan(L)
        L = svds(A,1,'largest','SubspaceDimension',60)^2 + 2*alpha;
    end

    % Group Lasso
    proj = @(v, tau) prox_group_lasso(v, tau, G);
    
    [ X_fista, iter, runtime ] = fista_general(grad, proj, Xinit, L, opts, calc_F)   ;

    h.iter = iter;
    h.runtime = runtime;

end


%% Proximal function

function x = prox_group_lasso(v, tau, groups)
    % v:        current point
    % tau:      lambda / L
    % groups:   cell array of index vectors, e.g. {G1, G2, ...}
    
    x = zeros(size(v));
    % Child
    for k = 1:numel(groups)
        G = groups{k};
        vg = v(G);
        ng = norm(vg, 2);
        s  = max(1 - tau/max(ng, eps), 0);   % handle ng = 0 safely
        x(G) = s * vg;
    end

    % y = x;
    % 
    % % Parent
    % groups = reshape(groups,10,55);
    % for i = 1:10
    %     G = sort( cell2mat(groups(i,:)) );
    %     vg = x(G);
    %     ng = norm(vg, 2);
    %     s  = max(1 - 0*tau/max(ng, eps), 0);   % handle ng = 0 safely
    %     x(G) = s * vg;
    % 
    % end
    % 
    % rel_err(x,y)

end