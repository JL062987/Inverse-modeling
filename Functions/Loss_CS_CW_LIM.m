function [Q,D] = Loss_CS_CW_LIM(Gam,A,B,N0,N1,dN0,At,Bt)

    [n,~,N] = size(A);
    dt = 1/N;
    tCp = (0:N-1)/N;
    tAp = tCp+0.5*dt;


    A1 = 1;
    A2 = 1;
    A3 = 1;
    A4 = 1e-5:0.1:1.5;


    % For NPP
    % A1 = -0.5:0.5:1.5;
    % A2 = -0.5:0.5:1.5;
    % A3 = -0.5:0.5:1.5;
    % A4 = +0.8:0.1:1.2;



    err_ModelSelection = nan(1,(length(A1)*length(A2)*length(A3)*length(A4)));
    flag_ModelSelection = nan(1,(length(A1)*length(A2)*length(A3)*length(A4)));
    Q_ModelSelection = cell(1,(length(A1)*length(A2)*length(A3)*length(A4)));
    D_ModelSelection = cell(1,(length(A1)*length(A2)*length(A3)*length(A4)));

    for i = 1:(length(A1)*length(A2)*length(A3)*length(A4))

        [a1,a2,a3,a4] = ReshapeIdx(i);

        weight = DefineWeight(A1(a1),A2(a2),A3(a3),n); % Define weight

        %%% Estimate noise covariance matrix for colored noise
        b = [];
        X = [];
        MM = Sym( N1 - pagemtimes(A,N0) );
        BB = B;
        for j = 1:N
            b = [ b ; diag(vec(weight)) * vec(MM(:,:,j)) ];
            X = [ X ; diag(vec(weight)) * ( kron(eye(n),BB(:,:,j)) + kron(BB(:,:,j),eye(n)) ) ];
        end

        Aieq = [];
        Bieq = [];
        Aeq = eye(n^2)-com_mat(n,n);
        beq = zeros(n^2,1);
        options = optimoptions('lsqlin','Display','off','MaxIterations',50);
        [ Q, ~, ~, exitflag ] = lsqlin(X,b,Aieq,Bieq,Aeq,beq,[],[],[],options);

        flag_ModelSelection(i) = double( exitflag < 0 );

        Q = reshape(Q,n,n);
        Q = Sym(Q)/2;
        Q = Q_Adjustment(Q,1);

        if norm(vec(Q)) < 1e-7 || sum(isnan(vec(Q))) > 0
            Q = 0.01*eye(n);
            disp("Q is artificially scaled.")
        end


        % The derivative of correlation function determined by A and Q
        N1new = pagemtimes(A,N0) + pagemtimes(Q,pagetranspose(B));


        % Estimate diffusion matrix for white noise by solving the FDR
        D = 0.5 * ( dN0 - Sym( N1new ) );
        D = Q_Adjustment(D,1); % Avoid negative eigenvalues

        % Estimate the model covariance
        Q_update = A4(a4)*Q;

        % The corresponding covariance function
        Dt = mat2fun(D,tAp);
        [Cxx,h_Period] = Periodic_Sol( @(t) kron( At(t), speye(n) ) + kron( speye(n), At(t) ) , @(t) vec( Sym( Q_update*Bt(t)' + Dt(t) ) ) , tCp );
        
        if h_Period.flag == 1 % Unable to solve covariance function
            z = nan;
        else
            Cxx = reshape(Cxx,n,n,[]);

            % Define weight
            w = DefineWeight(1,1,5,n);
            w = repmat( w, 1, 1, N );
            z = rel_err(w.*Cxx,w.*N0);

        end

        % Relative error between input and model covariance
        err_ModelSelection(i) = z;
        Q_ModelSelection{i} = Q_update;
        D_ModelSelection{i} = D;

        % A4(a4)
        % err_ModelSelection
    end

    [a,b] = min(err_ModelSelection);
    str = strcat( "FDR-residual: ", num2str(100*a,'%2.1f'), "%" );
    % disp( str )
    flag = flag_ModelSelection(b);
    Q = Q_ModelSelection{b};
    D = D_ModelSelection{b};

    % sq2QCetax = pagemtimes(Q,B2);
    % 
    % 
    % 
    % % The derivative of correlation determined by A and Q
    % N1new = pagemtimes(A,N0) + sq2QCetax;
    % % rel_err(N1,N1new)
    % 
    % % Estimate Q by solving the FDR
    % D = 0.5 * ( dN0 - Sym( N1new ) );
    % D = Q_Adjustment(D,1); % Avoid negative eigenvalues


    function [a1,a2,a3,a4] = ReshapeIdx(i)
        
        a1 = myquotient(i,length(A2)*length(A3)*length(A4));
        b1 = mymod(i,length(A2)*length(A3)*length(A4));
    
        a2 = myquotient(b1,length(A3)*length(A4));
        b2 = mymod(b1,length(A3)*length(A4));
    
        a3 = myquotient(b2,length(A4));
        a4 = mymod(b2,length(A4));
    end
            
end

function weight = DefineWeight(x,y,z,n) % Eq. (S10)

    weight = ones(n,n);
    weight(1,1:end) = x; 
    diagonalMask = logical(eye(n));
    weight(diagonalMask) = y; 
    weight(1,1) = z; 

    % Symmetrize
    weight = triu( weight );
    weight = weight + weight' - diag(diag(weight));
end


function y = Sym(x)
    y = x + pagetranspose(x);
end