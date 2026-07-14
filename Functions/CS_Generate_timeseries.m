function [xvec,h] = CS_Generate_timeseries(Type,A,Q,D,Gam,TT,dt,myopt)   

  
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
    

    h = struct();

    
    [n,~] = size(A(1)); % System dimension
    
    if Type == "White"
        % Q will be ignored.
        sq2D = @(t) sqrtm(2*D(t));
    else
        if isempty(D) || Type == "Colored"; D = []; end
        [A,~,sq2D] = CS_to_Aug_System(Gam,A,Q,D);
    end

    [ns,~] = size(A(1)); % Dimension of the (possible) augmented system


    Fsde = @(t,X) A(t)*X;
    Gsde = @(t,X) sq2D(t);

    T0 = 0; % Initial time
    X0 = randn(ns,1); % Initial condition
    T_spin = 10;
    nPeriods = (TT+T_spin)/dt; % Timespan

    obj = sde(Fsde,Gsde,'Correlation',eye(ns),'StartTime',T0,'StartState',X0);
    
    % Generate the time series
    [xvec,T] = simByEuler(obj,nPeriods,'DeltaTime',dt,'NTrials',1); % Simulation
   

    if Type == "White"
        xvec = xvec'; 
    else
        h.etavec = xvec(:,(n+1):end)'; % Store colored-noise process
        xvec = xvec(:,1:n)'; 
    end
    
    xvec = xvec(:,round(T_spin/dt+1):end-1); % Remove spin-up 

end