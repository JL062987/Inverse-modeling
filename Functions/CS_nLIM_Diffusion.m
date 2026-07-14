function Q = CS_Colored_nLIM_Diffusion(Type,B0,A0,C0,Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,Ktheta)

    [n,N] = size(E0);
    Dt = 1/N;


    Q = zeros(n,n,N);

    if Type == "White"

        for j = 1:N
        
            Bt = B0(:,:,:,j);
            At = A0(:,:,j);
            Ct = C0(:,j);
        
            M0t = M0(:,:,:,j);
            K0t = K0(:,:,j);
            E0t = E0(:,j);
        
            LHS = Ktheta(:,:,j);
        
            % Second-order
            [T1,x1] = Tensor_Mul(Bt,M0t,2);
            [T2,x2] = Tensor_Mul(At,K0t,1);
            [T4,x4] = Tensor_Mul(Ct,E0t,0);

            RHS = reshape( ( T1*x1 + T2*x2 + T4*x4 ), n, n );
            RHS = RHS + RHS';
            
            Q(:,:,j) = 0.5*( LHS - RHS );

            
        
        end

    else
        
        sqrt2QKetax = zeros(n,n,N);
        sqrt2QMetax = zeros(n,n,n,N);
        
        for j = 1:N
        
            Bt = B0(:,:,:,j);
            At = A0(:,:,j);
            Ct = C0(:,j);
        
            S0t = S0(:,:,:,:,j);
            M0t = M0(:,:,:,j);
            K0t = K0(:,:,j);
            E0t = E0(:,j);
        
            M1t = M1(:,:,:,j);
            K1t = K1(:,:,j);
        
            % Second-order
            [T1,x1] = Tensor_Mul(Bt,M0t,2);
            [T2,x2] = Tensor_Mul(At,K0t,1);
            [T4,x4] = Tensor_Mul(Ct,E0t,0);
        
            yK = vec( K1t ) - ( T1*x1 + T2*x2 + T4*x4 ); % sqrt2QKetax
            sqrt2QKetax(:,:,j) = reshape( yK, n, n );
        
            % Third-order
            [T11,x11] = Tensor_Mul(Bt,S0t,2);
            [T12,x12] = Tensor_Mul(At,M0t,1);
            [T14,x14] = Tensor_Mul(Ct,K0t,0);
        
            yM = vec( M1t ) - ( T11*x11 + T12*x12 + T14*x14 );
            sqrt2QMetax(:,:,:,j) = reshape( yM, n, n, n );
        
        
        end
        
        
        % Formulating Eq. (A8)
        
        for j = 1:N
            Bt = B0(:,:,:,j);
            At = A0(:,:,j);
        
            sqrt2QMetaxt = sqrt2QMetax(:,:,:,j);
            sqrt2QKetaxt = sqrt2QKetax(:,:,j);
        
            idx_f = mymod(j+1,N);
            idx_b = mymod(j-1,N);
        
            dLHSdt = ( sqrt2QKetax(:,:,idx_f) - sqrt2QKetax(:,:,idx_b) ) / (2*Dt) ;
            
            [T21,x21] = Tensor_Mul(sqrt2QMetaxt,permute(Bt,[2 3 1]),2);
            [T22,x22] = Tensor_Mul(sqrt2QKetaxt,At'-1/Gam*eye(n),1);
        
            Qt = reshape( Gam * ( vec( dLHSdt ) - (T21*x21+T22*x22) ), n, n );
            Qt = 0.5*(Qt+Qt');
        
            Q(:,:,j) = Qt;
        end
        
    
    end


    Q = mean(Q,3);
    
    Q = Q_Adjustment(Q,1);
    Q = 0.5*(Q+Q');
    

end