function Q = fourier_basis_real(M,n)

    if 2*n > M
        disp('n is too large')
        n = floor(M/2);
    end

    % Real Fourier basis over grid m=0..M-1 (one period).
    m = (0:M-1).';
    cols = ones(M,1)/sqrt(M);                  % DC
    for k = 1:n
        % cols = [cols, cos(2*pi*k*m/M)*sqrt(2/M)];
        % Nyquist: if M even and k==M/2, sin column is identically zero -> skip
        if (mod(M,2)==0 && k==M/2)
            cols = [cols, cos(2*pi*k*m/M)*sqrt(1/M)];
        else
            cols = [cols, cos(2*pi*k*m/M)*sqrt(2/M)];
            cols = [cols, sin(2*pi*k*m/M)*sqrt(2/M)];
        end
    end
    Q = cols;
    
end