function Q = CS_nLIM_Diffusion_Residual_Method(NoiseType,xvec,B0,A0,C0,Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta)

    [~,T] = size(xvec);
    [n,N] = size(E0);
    Dt = 1/N;
    
    
    DxDt = ( xvec(:,2:end) - xvec(:,1:end-1) ) / Dt ;
    
    Bx = zeros(n,T-1);
    Ax = zeros(n,T-1);
    Cx = zeros(n,T-1);
    for j = 1:(T-1)
        idx = mymod(j,N);
        Bx(:,j) = f(B0(:,:,:,idx),xvec(:,j));
        Ax(:,j) = A0(:,:,idx)*xvec(:,j);
        Cx(:,j) = vec(C0(:,idx));
    end
    
    % Qtrue = [1 0.5; 0.5, 1];
    % invsq2Q = inv( sqrtm(2*Qtrue) );
    
    Res = DxDt - (Bx+Ax+Cx);
    
    % Res = pagemtimes(invsq2Q,Res);
    % rel_err(Res,etavec(:,1:end-1))
    
    if NoiseType == "White"
    
        Q = Res*Res'/(T-1);
        Q = Q*Dt/2;
    
    else
        
        maxlag = 1;
        K = CorrelationFunction_ST_Obs(Res,maxlag);
        K = K(:,:,(maxlag+1):end);
        
        % size(K)
        Res_K0 = K(:,:,1);
        Res_K1 = ( K(:,:,2) - K(:,:,1) ) / Dt ;
        % logm( K(:,:,2) / K(:,:,1) ) / Dt
        
        
        [Q1, Q2] = ST_White_LIM(Res_K0,Res_K1);
        
        Q = -Gam*Q2/Q1;
    
    end
    
    Q = 0.5*(Q+Q');


end

%% 

function y = f(B,x)
    n = length(x);
    y = zeros(n,1);
    for i = 1:n
        for j = 1:n
            for k = 1:n
                y(i) = y(i) + B(i,j,k)*x(j)*x(k);
            end
        end
    end

end

function [A,Q] = ST_White_LIM(K,K1)
    A = K1/K;
    Q = -0.5*( A*K + K*A' );
end

function K = CorrelationFunction_ST_Obs(xvec,maxlag)

    % xvec: var * time (e.g. monthly data)
    % The lag variable is sampled as the same as input xvec.

    % Output
    % Stationary correlation function (lag-tau) covariance

    
    [n,~] = size(xvec);
    K = zeros(n,n,2*maxlag+1); % Observed correlation
    
    for i = 1:n
        for j = 1:n
            xi = xvec(i,:);
            xj = xvec(j,:);
            K(i,j,:) = xcorr(xi,xj,maxlag,"unbiased");
        end
    end

end