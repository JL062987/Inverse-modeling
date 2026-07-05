function [A,Q,D,h] = CS_LIM(Type,N0,N1,N2,Gam,myopt)

% Input
%   K: Cyclostationary correlation function
%   Gam: Noise correlation time

    if nargin < 6
        myopt = struct();
    end


    if isfield(myopt,'Solve_Diffusion')
        Solve_Diffusion = myopt.Solve_Diffusion;
    else
        Solve_Diffusion = "Loss_Function";
    end


    [n,~,N] = size(N0);
    dt = 1/N; % Timestep

    if isscalar(Gam); Gam = Gam*ones(n,1); end

    tCp = (0:(N-1))/N;
    tAp = tCp+0.5*dt;

    % Derivative of covariance
    dC = Phase_Derivative(N0) ;

    % Estimate the dynamics
    XX = zeros(n*N,n*N);
    ZZ = zeros(n,n*N);

    for i = 1:N
        idx = ((i-1)*n+1) : (i*n);

        % Diagonal
        XX(idx,idx) = Gam(1)*N1(:,:,i) + N0(:,:,i);

        % Phase derivatives
        XX(mymod(idx-n,n*N),idx) = -Gam(1)*N0(:,:,i)/(2*dt);
        XX(mymod(idx+n,n*N),idx) = Gam(1)*N0(:,:,i)/(2*dt);

        % LHS
        ZZ(:,idx) = Gam(1)*N2(:,:,i) + N1(:,:,i);
    end

    if isfield(myopt,'FM')
        AA = SolveEqwithFM(kron(XX',eye(n)),vec(ZZ),myopt.FM);
    else
        AA = reshape( kron(XX',eye(n)) \ vec(ZZ), n,n,[] );
    end

    A = reshape(AA,n,n,[]);
    
    % Check if the dynamcis is stable.
    if sum( abs(eig(Monodromy(A)))<1 ) ~= n % || CheckStability(A(1,1,:)) ||  CheckStability(A(n,n,:))

        h.flag = 1;
        
        Q = zeros(n,n);
        D = zeros(n,n,N);

    else

        h.flag = 0;

        if Type == "White"

            Q = zeros(n);

            D = 0.5*( dC - Sym( pagemtimes(A,myinterp1(tCp,N0,tAp)) ) );
            D = Q_Adjustment(D,1);
        
        else

            At = mat2fun(A,tAp);

            Gaminv = eye(n)./Gam;
            [BT,h_PerSol] = Periodic_Sol_demo( @(t) kron( -Gaminv + At(t), speye(n) ) , @(t) vec( Gaminv ) , tCp );
            BT = reshape(BT,n,n,[]); % This is B^T!!
            B = pagetranspose(BT); % Transpose back to be B!
            h.B = B;

            if h_PerSol.myflag == 1 || sum( isnan(B(:)) ) > 0 
                Q = zeros(n,n);
                D = zeros(n,n,N);
            else

                Bt = mat2fun(B,tCp); 
    
                %%% Estimate the noise covariance matrix Q

                if Solve_Diffusion == "Loss_Function"
                    if Type == "Colored"
                        [Q,~] = Loss_CS_Colored_LIM(Gam,A,B,N0,dC,At,Bt);
                        D = zeros(n,n,N);
                    elseif Type == "CW"    
                        [Q,D] = Loss_CS_CW_LIM(Gam,A,B,N0,N1,dC,At,Bt);    
                    end
                else

                    if Type == "Colored"
                        
                        MM = dC - Sym(pagemtimes( A, N0 )); % C' = AC + CA^T in Eq. (21), with the phase specifying in the third dimension
                        
                        b = zeros(n^2*N,1);
                        X = zeros(n^2*N,n^2);
                    
                        idx0 = 1:(n^2);
                        for j = 1:N
                            idx = idx0 + (j-1)*(n^2);
                            b(idx) = vec(MM(:,:,j));
                            X(idx,:) = ( kron(eye(n),B(:,:,j)) + kron(B(:,:,j),eye(n)) );
                        end

                        Q = SolveEqwithSymmetric(X,b);


                        D = zeros(n,n,N);

                    elseif Type == "CW"
                        
    
                        % Diffusion matrix for colored noise
                        b = [];
                        X = [];
                        
                        MM = Sym( N1 - pagemtimes(A,N0) );

                        for j = 1:N
                            b = [ b ; vec(MM(:,:,j)) ];
                            X = [ X ; ( kron(eye(n),B(:,:,j)) + kron(B(:,:,j),eye(n)) ) ];
                        end
           
                        Q = SolveEqwithSymmetric(X,b);

                       
                        % Diffusion matrix for white noise
                        D = 0.5*( dC - Sym( pagemtimes(A,N0) + pagemtimes(Q,pagetranspose(B)) ) );
                        % D = 0.5*( dC-Sym(N1) );
                        
                    
                         
                    end
                end
            end
        end
    end
    
    
    function x = SolveEqwithFM(A,b,k)

        if nargin == 2
            k = 0;
        end

        % Solve Ax = b where x is subjected to Fourier truncation
        m = length(b)/N; % Number of variables
    
        % Rearrangement
        p = vec(reshape(1:(m*N),m,[])'); % Permutation vector
        [~,q] = sort(p); % Inverse permutation
    
        Ap = A(p,p);
        bp = b(p);
    
        Q0 = fourier_basis_real(N,k);
        F = myBlkdiag(Q0,m);
    
        % A = A*F;
        x = F*((Ap*F)\bp);
        x = x(q);
    end

    function Q = SolveEqwithSymmetric(X,b)

        Aieq = [];
        Bieq = [];
        Aeq = eye(n^2)-com_mat(n,n); % Make sure that Q is symmetric when solved by lsqlin 
        beq = zeros(n^2,1);
        options = optimoptions('lsqlin','Display','off','MaxIterations',100);
        Q = lsqlin(X,b,Aieq,Bieq,Aeq,beq,[],[],[],options);
        Q = reshape(Q,n,n); 

        % In case of numerical errors
        Q = Sym(Q)/2;
        Q = Q_Adjustment(Q,1);
    end



end


function y = Sym(x)
    y = x + pagetranspose(x);
end


function dCC = Phase_Derivative(CC)

    % A: the time dependent dynamics
    % C: the covariacne function
    % dCC: the time derivative of C
    
    
    [~,~,T] = size(CC);
    
    dt = 1/T; 
    
    dCC = zeros(size(CC));

    for t = 1:T  

        t1 = mymod( t, T );
        t2 = mymod( t+1, T ); 
        dCC(:,:,t) = squeeze( CC(:,:,t2) - CC(:,:,t1) ) / (dt); 

    end
    
end

function z = CheckStability(X)
% Check if the high frequency part
    N = length(X);
    Y = fft(vec(X));
    P2 = abs(Y/N);
    P1 = P2(1:N/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    z = P1(end) > P1(1); % z = 1: Unstable
end
