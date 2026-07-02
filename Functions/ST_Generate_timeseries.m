function [xvec,h] = ST_Generate_timeseries(Type,A,Q,D,Gam,TT,dt,myopt)   

  
    % A: Dynamical matrix
    % D: Diffusion matrix
    % Gam: Noise correlation time
    % TT: Desired time span
    % dt: Sampling interval
    
    %%% The synthetic data
    % xvec: State variables
    % etavec: Colored noise process
    

    if nargin < 8
        myopt = struct();
    end
    
    if isfield(myopt,"Method")
        Method = myopt.Method;
    else
        Method = "Euler";
    end

    h = struct();



    
    [n,~] = size(A); % System dimension
    
    if Type == "White"
        % Q will be ignored.
        sq2D = sqrtm(2*D);
    else
        if isempty(D) || Type == "Colored"; D = zeros(n); end
        [A,~,sq2D] = ST_to_Aug_System(Gam,A,Q,D);
    end
    
    
    [ns,~] = size(A); % Dimension of the augmented system (= n for white; = 2n for colored or CW)
    
    DriftRate = drift(zeros(ns,1),A); % Diffusion rate 
    DiffusionRate = diffusion(zeros(ns,1),sq2D); % Diffusion rate
    T0 = 0; % Initial time
    X0 = zeros(ns,1); % Initial condition
    SDEDDO = sdeddo(DriftRate,DiffusionRate,'Correlation',eye(ns),'StartTime',T0,'StartState',X0); % Sde model
    nPeriods = (TT+10)/dt; % Timespan

    if Type == "Colored" && Method == "Milstein2"
        [xvec,~] = simByMilstein2(SDEDDO,nPeriods,'DeltaTime',dt,'NTrials',1); % Simulation
    else % simByMilstein2 requires diagonal diffusion matrix
        [xvec,~] = simByEuler(SDEDDO,nPeriods,'DeltaTime',dt,'NTrials',1); % Simulation
    end

    if Type == "White"
        xvec = xvec'; 
    else
        h.etavec = xvec(:,(n+1):end)'; % Store colored-noise process
        xvec = xvec(:,1:n)'; 
    end
    
    xvec = xvec(:,round(10/dt+1):1:end-1);

end