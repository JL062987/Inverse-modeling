function [z0,B,A,C,h] = CS_nLIM_Dynamics(Type,Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,idx_nonzeroB,idx_nonzeroA,idx_nonzeroC,myopts)

    %%% Input
    % Gam: Noise correlation time
    % Statistics for x.
    % idx_nonzero: Specify nonzero (active) components in the dynamics
    % myopts (optional): Specify the use of Fourier mode truncation

    %%% Output
    % z0: The estimated parameters stacked into a column vector
    % B, A, C: The estimated parameters in tensor forms


    if nargin < 16
        myopts = struct();
        myopts.Method = "Fixed_phase";
    end

    if Type == "White"; Gam = 0; end

    [n,N] = size(Etheta);

    % Prelocate variables
    B = zeros(n,n,n,N);
    A = zeros(n,n,N);
    C = zeros(n,N);

    % The active components (those not prescribed to be zero)
    idx_nonzeroB = logical(idx_nonzeroB);
    idx_nonzeroA = logical(idx_nonzeroA);
    idx_nonzeroC = logical(idx_nonzeroC);
    nB = sum(idx_nonzeroB);
    nA = sum(idx_nonzeroA);
    nC = sum(idx_nonzeroC);
    
    idx_nonzero = logical( [ idx_nonzeroB; idx_nonzeroA; idx_nonzeroC ] );
    Nv_nonzero = sum(idx_nonzero);
   


    % With T0 and b0, the dynamical parameters ("V" in Eq. (8)) are solved
    % by matrix inverse using lsqr, an iterative matlab built-in function
    tol = 1e-6; % Tolorence
    maxit = 10000; % Max iteration
    
    if myopts.Method == "Fixed_phase"

        [Tp1,bp1,iq] = CS_nLIM_Formulate_Eq_6(Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,idx_nonzero);
        z0 = lsqr(Tp1,bp1,tol,maxit); % Arrange as [B(1),...,B(end),A(1),...,A(end),C(1),...,C(end)]
        z0 = z0(iq); % Arrange as [B(1),A(1),C(1),...,B(end),A(end),C(end)]


    elseif myopts.Method == "Phase_average"

        [Tp,bp,iq] = CS_nLIM_Formulate_Eq_6(Gam,E0,K0,K1,K2,M0,M1,M2,S0,S1,Etheta,idx_nonzero);
        
        Q0 = fourier_basis_real(N,myopts.FM); % Projection matrix
        Q1 = myBlkdiag(Q0,nB+nA+nC);

        [z0,~] = lsqr(Tp*Q1,bp,tol,maxit);
        z0 = Q1*z0; % Projection back
        z0 = z0(iq); 

    else
        disp('Wrong myopts.Method!')
    end


    % Transform the dynamical parameters from a single long vector back into the tensor form
    idx0 = reshape( 1:(Nv_nonzero*N), Nv_nonzero, [] );
    for j = 1:N
        idx   = idx0(:,j);

        B0 = zeros(n^3,1);
        B0(idx_nonzeroB) = z0(idx(1:nB));
        B(:,:,:,j) = reshape(B0,n,n,n);

        A0 = zeros(n^2,1);
        A0(idx_nonzeroA) = z0(idx((nB+1):(nB+nA)));
        A(:,:,j) = reshape(A0,n,n);

        C0 = zeros(n,1);
        C0(idx_nonzeroC) = z0(idx((nB+nA+1):(nB+nA+nC)));
        C(:,j) = C0;

    end

    h = struct();

end