function [xvec,h] = CS_Nonlinear_Generate_timeseries(Type,F,Q,Gam,TT,dt,opt)

    % Do not support CW type!

    if nargin < 7
        opt = struct();
    end

    if ~isfield(opt,'NoiseGeneration')
        opt.NoiseGeneration = "augmented";
    end


    T_spin = 10;

    [n,~] = size(Q); % Dimension of the (possible) augmented system

    sq2Q = sqrtm(2*Q);

    if Type == "Colored" && opt.NoiseGeneration == "precompute"

        if isfield(opt,'eta')
            eta = opt.eta;
        else
            eta = zeros(n,(TT+T_spin)/dt+1);
            for j = 1:n
                eta(j,:) = Generate_Colored_Noise(Gam,TT+T_spin,dt); 
            end

        end
        

        x0 = zeros(n,1);

        if ~round((TT+T_spin)/dt+1) == length(eta)
            disp('Do you include spin-up period [0,10]?')
        end

        xvec = SDE_integrate_color(F,Q,T_spin+TT,x0,dt,eta);

        h = struct();
    
    else

        if Type == "White"
    
            Fsde = F;
            if isa(Q,'function_handle')
                Gsde = @(t,X) sqrtm(2*Q(t));
            else
                sq2Q = sqrtm(2*Q);
                Gsde = @(t,X) sq2Q;
            end
            ns = n;

        else
    
            Gaminv = 1/Gam;
            sq2Q_aug = [zeros(n), zeros(n); zeros(n), Gaminv*eye(n)];
        
            Fsde = @(t,X) [ F(t,X(1:n)) + sq2Q*X(n+1:end); -Gaminv*X(n+1:end) ];
            Gsde = @(t,X) sq2Q_aug;
    
            ns = 2*n;

        end

        T0 = 0; % Initial time
        X0 = randn(ns,1); % Initial condition
        nPeriods = (TT+T_spin)/dt; % Timespan
    
        obj = sde(Fsde,Gsde,'Correlation',eye(ns),'StartTime',T0,'StartState',X0);
        
        % Generate the time series
        [xvec,T] = simByEuler(obj,nPeriods,'DeltaTime',dt,'NTrials',1); % Simulation
        xvec = xvec'; 
    end
   

    if Type == "White"
        % Do nothing
    else
        % h.etavec = xvec(:,(n+1):end)'; % Store colored-noise process
        xvec = xvec(1:n,:); 
    end
    
    i0 = round(10/dt) + 1; % Remove spin-up 
    xvec = xvec(:,i0:end-1);
    
end
