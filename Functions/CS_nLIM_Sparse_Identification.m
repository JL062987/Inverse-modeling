function [idx_nonzeroB,h] = CS_nLIM_Sparse_Identification(Type,Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,opts_in)

    if isfield(opts_in,"Method")
        Method = opts_in.Method;
    else
        Method = "Fixed_phase";
        disp('Default method: Fixed_phase.')
    end
    
    if ~isfield(opts_in,"Solve_Lambda")
        opts_in.Solve_Lambda = "Grid_Search";
        opts_in.user_opts = struct();
    end

    if ~isfield(opts_in,"FISTA_opt")
        FISTA_opts = struct();
    else
        FISTA_opts = opts_in.FISTA_opts;
    end

    if Type == "White"; Gam = 0; end

    %%%

    [n,N] = size(E0);


    %%% Formulate Eq. (6)
    [II,JJ,KK] = ndgrid(1:n,1:n,1:n);
    idx = [ vec( JJ <= KK ); true(n^2,1); true(n,1) ];

    [Tp,bp] = CS_nLIM_Formulate_Eq_6(Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,idx);

    %%% Formulate Eq. (7)
    
    neff = n^2*(n+1)/2;
    [Atilde_lasso,btilde_lasso,~] = Schur_Complement(Tp,bp,neff*N); % 

    %%% Apply Fourier-mode truncation or not
    if Method == "Phase_average"

        if isfield(opts_in,"FM")
            FM = opts_in.FM; % The number of Fourier modes
        else 
            FM = 2;
        end
        Q0 = fourier_basis_real(N,FM); % Projection matrix
        FMnum = size(Q0,2); % Const + sum( sin + cos )
        Q1 = myBlkdiag(Q0,neff);
        Atilde_lasso = Atilde_lasso*Q1; 

    else % Fixed-phase
        
        FMnum = N;

    end

    % isfield(inputopts,"lambda")
    
    
    %%% Objective function for lasso regression
    idx_lasso = reshape(1:(neff*FMnum),[],neff);  
    S = cell(neff,1);
    for i = 1:neff; S{i} = idx_lasso(:,i); end % Lasso group
    

    FISTA_opts.tol = 1e-5; 
    ObjFcn = @(x) StandardizeLasso(Atilde_lasso,btilde_lasso,S,x,FMnum,neff,FISTA_opts); 


    %%% Set up parameters
    if opts_in.Solve_Lambda ~= "Grid_Search"
        Scheme = @opts_in.Solve_Lambda;
        user_opt = opts_in.user_opts;

    else
        Scheme = @Grid_search_lambda;
        user_opt = opts_in.user_opts;

        if ~isfield(user_opt,"lambda_range")
            if Method == "Phase_average"
                user_opt.lambda_range = 0:1:30;
            else
                user_opt.lambda_range = 0:0.1:5;
            end
        end
        
    end

    % Find out the best Lambda
    [idx,h] = Scheme(ObjFcn,user_opt);
    idx_nonzeroB = zeros(n,n,n);
    idx_nonzeroB(JJ <= KK) = idx; % From vectorized and non-repeated index to the standard (i,j,k) index

    logical_idx = logical(idx_nonzeroB);
    h.ActiveComponents = [vec(II(logical_idx)) vec(JJ(logical_idx)) vec(KK(logical_idx))];
    

end

%% 


function [idx,h] = Grid_search_lambda(f,user_opt)

    lambda_range = user_opt.lambda_range; % Range for the grid search

    iter_count = zeros(1,length(lambda_range));
    t_count    = zeros(1,length(lambda_range));
    Lasso_path = cell(1,length(lambda_range));

    for j = 1:length(lambda_range)
        lambda = lambda_range(j);
        [x1,h] = f(lambda);
        iter_count(j) = h.iter; % Iteration count in FISTA
        t_count(j) = h.runtime; % Runtime 
        Lasso_path{j} = vec( vecnorm(x1) ); % Lasso path
    end

    Lasso_path = cell2mat(Lasso_path);
    Lasso_path( Lasso_path < 1e-2 ) = 0; % Mask out zero 
    [~,lambda_idx] = min( sum( Lasso_path ~= 0, 1 ) );
    
    idx = Lasso_path(:,lambda_idx)>0; % Indicator of the active components in B.

    % Store results
    h.Lambda_range = lambda_range;
    h.Lasso_path = Lasso_path;
    h.iter_count = iter_count;
    h.t_count = t_count;

end

function [ x, h ] = StandardizeLasso(A,b,S,lambda,FMnum,neff,opts)

    % 1. Center b
    b_mean = mean(b);
    b_standardize = b - b_mean;
    
    % 2. Center and scale A
    A_mean = mean(A, 1);              % mean of each column (1 x n)
    A_centered = A - A_mean;          % subtract mean
    A_std = std(A_centered, 0, 1);    % standard deviation of each column (1 x n)
    A_standardize = A_centered ./ A_std;     % elementwise division
    
    [beta,h_FISTA] = CS_nLIM_Group_Lasso_FISTA(sparse(A_standardize),b_standardize,S,lambda,opts);
    
    % 3. Unstandardize the coefficients
    x = beta ./ A_std';  % (n x 1)
    % intercept = b_mean - A_mean * x;  % recover intercept

    x = reshape(x,FMnum,neff);

    h.std = A_std;

    h.iter = h_FISTA.iter;
    h.runtime = h_FISTA.runtime;
    
end


