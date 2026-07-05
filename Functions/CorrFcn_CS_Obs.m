function K = CorrFcn_CS_Obs(x,dt,maxlag)

    % x: State vectors
    % dt: sampling timestep
    % maxlag: maxlag for the ouput lagged covariance function

    [n,m] = size(x);
    N = round(1/dt);
    P = m/N;   % number of periods
    
    % [var,phase,period]
    X = reshape(x,n,N,[]);
   
    K = zeros(n,n,N,maxlag+1); 

    for t = 1:N % Run over phase
        for lag = 0:maxlag

            tlag = round(mymod(t+lag,N)); % phase of t + lag
            plag = round(myquotient(t+lag,N));

            L = P - plag + 1; % Number of usable samples

            X_future = reshape(X(:,tlag,plag:P),n,L);
            X_now    = reshape(X(:,t,1:L),n,L);

            K(:,:,t,lag+1) = X_future * X_now' / L;

        end
    end

end