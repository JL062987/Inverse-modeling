function varargout = CorrFcn_CS_Obs(x,dt,maxlag,computeHigherOrder)

    % x: State vectors
    % dt: sampling timestep
    % maxlag: maxlag for the ouput lagged covariance function

    if nargin < 4
        computeHigherOrder = false;
    end

    [n,T] = size(x);
    N = round(1/dt);
    P = T/N;   % number of periods


    if mod(T,N) ~= 0
        error('The time-series length must be divisible by N = 1/dt.');
    end
    
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

    if computeHigherOrder

    
        X = reshape(x,n,N,P);
    
        E = mean(X,3);
    
        M = zeros(n,n,n,N,maxlag+1,'like',x);
        S = zeros(n,n,n,n,N,maxlag+1,'like',x);
    
        for t = 1:N
    
            % Data at reference phase
            Xt = reshape(X(:,t,:),n,P);
    
            % Calculate these only once for this phase
            P2 = reshape( ...
                reshape(Xt,n,1,P) .* ...
                reshape(Xt,1,n,P), ...
                n^2,P);
    
            P3 = reshape( ...
                reshape(Xt,n,1,1,P) .* ...
                reshape(Xt,1,n,1,P) .* ...
                reshape(Xt,1,1,n,P), ...
                n^3,P);
    
            for lag = 0:maxlag
    
                targetPhase = mod(t-1+lag,N) + 1;
                periodShift = floor((t-1+lag)/N);
    
                nSample = P-periodShift;
    
                if nSample <= 0
                    M(:,:,:,t,lag+1) = NaN;
                    S(:,:,:,:,t,lag+1) = NaN;
                    continue
                end
    
                Xlag = reshape( ...
                    X(:,targetPhase,periodShift+1:P), ...
                    n,nSample);
    
                M(:,:,:,t,lag+1) = reshape( ...
                    Xlag * P2(:,1:nSample).' / nSample, ...
                    n,n,n);
    
                S(:,:,:,:,t,lag+1) = reshape( ...
                    Xlag * P3(:,1:nSample).' / nSample, ...
                    n,n,n,n);
            end
        end
    end

    % Select outputs

    if computeHigherOrder

        if nargout ~= 4
            error(['When computeHigherOrder is true, use four outputs: ',...
                   '[E,K,M,S] = CorrFcn_CS_Obs(...)']);
        end


        varargout = {E,K,M,S};

    else

        if nargout ~= 1
            error(['When computeHigherOrder is false, use one output: ',...
                   'K = CorrFcn_CS_Obs(...)']);
        end

        varargout = {K};
    end

end