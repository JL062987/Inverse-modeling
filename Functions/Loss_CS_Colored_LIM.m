function [Q,h] = Loss_CS_Colored_LIM(Gam,A,B,N0,dN0,At,Bt)



    [n,~,N] = size(A);
    tCp = (0:N-1)/N;
    dt = 1/N;

    % Weight for the Frobenius norm
    weight0 = DefineWeight(1,3,10,n);
    
    A1 = 1:1:1;
    A2 = 1:1:1;
    A3 = 1:1:5;

    A4 = 0.7:0.05:1.0;

    err_ModelSelection = nan(1,(length(A1)*length(A2)*length(A3)*length(A4)));
    Q_ModelSelection = zeros(n,n,(length(A1)*length(A2)*length(A3)*length(A4)));

    for i = 1:(length(A1)*length(A2)*length(A3)*length(A4))
    
        [a1,a2,a3,a4] = ReshapeIdx(i);
        
        % Define weight as Eq. (S10)
        weight = DefineWeight(A1(a1),A2(a2),A3(a3),n);
  
        % Solve Q via Eq. (S11)
        MM = dN0 - Sym(pagemtimes( A, N0 )); % C' = AC + CA^T in Eq. (21), with the phase specifying in the third dimension

        b = zeros(n^2*N,1);
        X = zeros(n^2*N,n^2);

        idx0 = 1:(n^2);
        for j = 1:N
            idx = idx0 + (j-1)*(n^2);
            b(idx) = diag(vec(weight)) * vec(MM(:,:,j));
            X(idx,:) = diag(vec(weight)) * ( kron(eye(n),B(:,:,j)) + kron(B(:,:,j),eye(n)) );
        end

        Aieq = [];
        Bieq = [];
        Aeq = eye(n^2)-com_mat(n,n); % Make sure that Q is symmetric when solved by lsqlin 
        beq = zeros(n^2,1);
        options = optimoptions('lsqlin','Display','off');
        [ Q3, ~, ~, exitflag ] = lsqlin(X,b,Aieq,Bieq,Aeq,beq,[],[],[],options);
        if exitflag ~= 1
            disp("Check the potential numerical issue when using lsqlin!!")
        end
        Q3 = reshape(Q3,n,n);
        
        % Rescale Q, as mentioned in Line 85 in the supplementary text
        Q = A4(a4)*Q_Adjustment(Q3,1);

        % Reproduced covariance
        [C,h_Period] = Periodic_Sol( @(t) kron( At(t), speye(n) ) + kron( speye(n), At(t) ) , @(t) vec( Sym(Q*Bt(t)') ) , tCp );

        if h_Period.flag == 1 % Unable to solve covariance function
            z = nan;
        else
            C = reshape(C,n,n,[]);
        
            % Relative error
            w = repmat( weight0, 1, 1, N );
            z = rel_err(w.*C,w.*N0);
        end

        err_ModelSelection(i) = z;
        Q_ModelSelection(:,:,i) = Q;

    end

    % Pick up minimizer Q corresponding to the minimizer (a1,a2,a3,a4) shown below Line 86 in the supplementary text
    [~,b] = min(err_ModelSelection);
    Q = Q_ModelSelection(:,:,b);
   
    % The minimizer (a1,a2,a3,a4)
    h = struct();
    [h.a1,h.a2,h.a3,h.a4] = ReshapeIdx(b);




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