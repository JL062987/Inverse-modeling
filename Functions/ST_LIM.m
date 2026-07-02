function [A,Q,D] = ST_LIM(Type,K0,K1,K2,Gam)

    % Type: specify the noise type
    % K0, K1, K2: lag-0 covariance, first-order and second-order derivatives of the lagged covariance
    % Gam: prescribed noise correlation time (a constant)

    % A: Dynamical matrix for white noise
    % Q: Stochastic matrix for colored noise
    
    Sym = @(x) x+x';
    [n,~] = size(K0); 

    if Type == "White"; Gam = 0; end
    
    A = ( Gam*K2 + K1 ) / ( Gam*K1 + K0 );
    
    if Type == "White"
    
        Q = 0;
        D = -0.5*Sym(A*K0);
    
    elseif Type == "Colored" || Type == "CW"
    
        In = eye(n);
        B = inv( In - Gam*A );
        Vec2Mat = @(x) reshape(x,n,n); % Convert a vectorized matrix into its original form
    
        if Type == "Colored"
            Q = Vec2Mat( -(kron(B,In)+kron(In,B)) \ vec(Sym(A*K0)) ); 
            D = 0;
        else % CW
            Q = Vec2Mat( -(kron(B,In)+kron(In,B)) \ vec(Sym(A*K0-K1)) );
            Q = Q_Adjustment(Q,1); % Eigenvalue adjustment
    
            K1new = A*K0 + Q*B'; % The first-order derivative of lagged covariance determined by LIM

            D = -0.5*Sym(K1new);
            D = Q_Adjustment(D,1); % Eigenvalue adjustment

        end
        
    end


end