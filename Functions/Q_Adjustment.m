function Qnew = Q_Adjustment(Q,Adjustment_type)

    if nargin < 2
        Adjustment_type = 1;
    end

    [~,~,N] = size(Q);
    Qnew = zeros(size(Q));

    for j = 1:N
        [U,D] = eig( squeeze( Q(:,:,j) ) );
        D = real(diag(D));
        D0 = D;
        D0(D0<0) = 1e-7; % Remove negative eigenvalues
    
        % Scale eigenvalues
        switch Adjustment_type
            case 1
                D = diag( abs(sum(D)/sum(D0)) * D0 ); % Trace-preserving
            case 2
                D = diag( norm(D,1)/norm(D0,1) * D0 );
            case 3
                D = diag( norm(D,2)/norm(D0,2) * D0 );
        end

        % Reconstruct Q
        Qnew(:,:,j) = real( U*(D /U) );

    end

end