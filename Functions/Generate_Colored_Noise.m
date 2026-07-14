function eta = Generate_Colored_Noise(Gam,TT,dt)

    % Gam: Noise correlation time
    % TT: Total timespan
    % dt: Timestep

    T0 = 0; % Start time
    X0 = 0; % Initial condition

    Gaminv = 1/Gam;

    Fsde = @(t,X) -Gaminv*X;
    Gsde = @(t,X) Gaminv;
    SDE = sde(Fsde,Gsde,'Correlation',1,'StartTime',T0,'StartState',X0);
    
    % The number of query points 
    nPeriods = (TT+10)/dt; % +10 spin-up time
    
    % Generate the time series
    [eta,~] = simByEuler(SDE,nPeriods,'DeltaTime',dt,'NTrials',1);
    
    idx = (10/dt+1):(nPeriods+1); % Remove the span-up time
    eta = eta(idx)'; % Pick up the colored noise forcing

end