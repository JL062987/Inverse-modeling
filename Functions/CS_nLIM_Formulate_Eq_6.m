function [Tp,bp,q] = CS_nLIM_Formulate_Eq6(Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,idx_nonzero)


    [n,N] = size(E0);
    Dt = 1/N;

    idx_nonzero = vec( repmat( idx_nonzero, N, 1) );

    Nv = n^1+n^2+n^3; % Number of variables in a snapshot
    
    T0 = zeros( Nv*N, Nv*N );
    b0 = zeros( Nv*N, 1);
    
    idx0 = reshape( 1:(Nv*N), Nv, [] );
    
    for j = 1:N
    
        % Some zero tensors as fillers
        O3 = zeros(n,n,n);
        O2 = zeros(n,n);
        O1 = zeros(n,1);
        
        %%% Second-order
        % Terms irrelevant to the derivatives of the dynamical parameters
        SymM1 = M1(:,:,:,j)+permute(M1(:,:,:,j),[2 1 3]);
        [TB2,~] = Tensor_Mul( O3, Gam*SymM1 + M0(:,:,:,j), 2 );
        [TA2,~] = Tensor_Mul( O2, Gam*K1(:,:,j) + K0(:,:,j), 1 );
        [TC2,~] = Tensor_Mul( O1, E0(:,j), 0 );
    
        % Terms involving to the derivatives of the dynamical parameters
        [T7,~] = Tensor_Mul( O3,M0(:,:,:,j),2);
        [T8,~] = Tensor_Mul( O2,K0(:,:,j),1);
        [T9,~] = Tensor_Mul( O1,E0(:,j),0);
    
        b2 = vec( Gam*K2(:,:,j) + K1(:,:,j) ) ;
    
        %%% Third-order
        % Terms irrelevant to the derivatives of the dynamical parameters
        SymS1 = S1(:,:,:,:,j)+permute(S1(:,:,:,:,j),[2 1 3 4]);
        [TB3,~] = Tensor_Mul(O3,Gam*SymS1+S0(:,:,:,:,j),2);
        [TA3,~] = Tensor_Mul(O2,Gam*M1(:,:,:,j)+M0(:,:,:,j),1);
        [TC3,~] = Tensor_Mul(O1,K0(:,:,j),0);
    
        % Terms involving to the derivatives of the dynamical parameters
        [T17,~] = Tensor_Mul(O3,S0(:,:,:,:,j),2);
        [T18,~] = Tensor_Mul(O2,M0(:,:,:,j),1);
        [T19,~] = Tensor_Mul(O1,K0(:,:,j),0);
    
        b3 = vec( Gam*M2(:,:,:,j) + M1(:,:,:,j) ) ;
    
    
        % First-order
        [TB1,~] = Tensor_Mul(O3,K0(:,:,j),2);
        [TA1,~] = Tensor_Mul(O2,E0(:,j),1);
        TC1 = eye(n);
        
        b1 = vec(Etheta(:,j));    
    
        % Augmented
        Tt = [ TB3 TA3 TC3 ; TB2 TA2 TC2 ; TB1 TA1 TC1 ];
        bt = [ b3 ; b2 ; b1 ];
    
        % Overall
        idx_b = idx0(:,mymod(j-1,N));
        idx   = idx0(:,j);
        idx_f = idx0(:,mymod(j+1,N));
        
        Tt_corr = Gam*[ T17 T18 T19 ; T7 T8 T9 ; zeros(n,Nv) ]/(2*Dt);
        
        T0(idx,idx)   = Tt;
        T0(idx,idx_b) = -Tt_corr;
        T0(idx,idx_f) = +Tt_corr;
        b0(idx) = bt;
    
    end

    
    % Remove redundant linear equations
    T = T0(idx_nonzero,idx_nonzero); b = b0(idx_nonzero);    
    
    % Permute such that the state variables aligns as [ B(1) ... B(end); A(1),...A(end) ; C(1)...C(end) ].
    p = vec( reshape( 1:length(b), [], N )' ); % Permutation
    [~,q] = sort(p); % Inverse permutation
    
    Tp = T(p,p); bp = b(p); % The matrix and vector of Eq. (6)
    

end