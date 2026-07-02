function [K,h] = ST_CorrFcn(Type,dt,maxlag,A,Q,D,Gam)

    % Type: specify noise type
    % dt: timestep
    % maxlag: max time-lag considered
    % A: linear dynamics
    % Q: diffusion matrix for colored noise
    % D: diffusion matrix for white noise
    % Gam: noise correlation time (a constant)

    % Output
    % K: Lagged covariance function sampled at (0:maxlag)*dt
    % h: History
    % ---- Kaug: The lagged covariance function for the augmented system for colored or CW case
    
    if Type == "White"

        if nargin < 6
            D = Q; % If D is missing, then treat Q as the diffusion matrix for white noise.
        end

        K = fun(A,D,maxlag,dt);
        h = struct();
    
    elseif Type == "Colored" || Type == "CW"

        [n,~] = size(A);
        if isempty(D) || Type == "Colored"; D = zeros(n); end 
        if isscalar(Gam); Gam = Gam*ones(n,1); end

        [L,S] = ST_to_Aug_System(Gam,A,Q,D); % Augmented dynamical and diffusion matrices

        % L = [ A, sqrtm(2*Q); zeros(n), -diag(1./Gam) ]; % Augmented dynamical matrix
        % S = [ D, zeros(n);   zeros(n), diag( 1./(2*Gam.^2) )]; % Augmented stochastic matrix
        
        K = fun(L,S,maxlag,dt);
        h.Kaug = K;
        K = K(1:n,1:n,:);
        
    end


end

function K = fun(Dyn,Sto,maxlag,dt)

    [n,~] = size(Dyn);
    Vec2Mat = @(x) reshape(x,n,n);
    In = eye(n);
    tau = (0:maxlag)*dt; % Lag variable
    K = zeros(n,n,length(tau));
    
    % Computation
    K(:,:,1) = Vec2Mat( -(kron(Dyn,In)+kron(In,Dyn)) \ vec(2*Sto) ); % Covariance
    expA = expm(Dyn*dt);
    for j = 2:length(tau)
        K(:,:,j) = expA * K(:,:,j-1);
    end


end