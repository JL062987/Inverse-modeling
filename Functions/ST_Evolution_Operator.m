function [G,Gd,Gs] = ST_Evolution_Operator(Type,C,A,Q,D,Gam,tlead)

    % D does not involve in the evolution operator!

    [n,~] = size(A);
        
    Gd = zeros(n,n,length(tlead));
    Gs = zeros(n,n,length(tlead));
    
    Vec2Mat = @(x) reshape(x,n,n);
    
    if Type == "Colored" || Type == "CW"
        if isscalar(Gam); Gam = Gam*ones(n,1); end
        B = Vec2Mat( ( eye(n^2) - kron( diag(Gam), A ) ) \ vec(eye(n)) )  ;
        sqrt2Q = sqrtm(2*Q);
        Coupling = (sqrt2Q/2)*B'; % Coupling
        L = [A, sqrt2Q; zeros(n), diag(-1./Gam)]; % The dynamical matrix for the augmented system
        sd = 2*n; % System dimension
    elseif Type == "White"
        L = A;
        sd = n;
    end
        
        
    if sum(abs(diff(diff(tlead)))) > 1e-7 || length(tlead) < 3
        for j = 1:length(tlead)
            t = tlead(j);
            Evo_Matrix = expm(L*t);
            if Type == "Colored" || Type == "CW"
                Gd(:,:,j) = Evo_Matrix(1:n,1:n); % Deterministic contribution
                Gs(:,:,j) = Evo_Matrix(1:n,(n+1):(2*n))*(Coupling/C); % Stochastic contribution
            elseif Type == "White"
                Gd(:,:,j) = Evo_Matrix;
            end
        end
    else
        Evo_Matrix = eye(sd);
        dt = tlead(2) - tlead(1);
        M = expm(L*dt);
        for j = 1:length(tlead)

            Evo_Matrix = M*Evo_Matrix;

            if Type == "Colored" || Type == "CW"
                Gd(:,:,j) = Evo_Matrix(1:n,1:n); % Deterministic contribution
                Gs(:,:,j) = Evo_Matrix(1:n,(n+1):(2*n))*(Coupling/C); % Stochastic contribution
            elseif Type == "White"
                Gd(:,:,j) = Evo_Matrix;
            end
        end
    end

    G = Gd + Gs; % Total


end