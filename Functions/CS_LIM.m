function [A,Q,D,h] = CS_LIM(Type,N0,N1,N2,Gam,myopt)

    % Input
    %   K: Cyclostationary correlation function
    %   Gam: Noise correlation time

    % Output
    %   A, Q, D: System parameters
    %   h (struct): history

    h = struct();

    if nargin < 6
        myopt = struct();
    end


    if isfield(myopt,'Solve_Diffusion')
        Solve_Diffusion = myopt.Solve_Diffusion;
    else
        Solve_Diffusion = "Loss_Function";
    end

    if ~isfield(myopt,'Fourier_mode_when_unstable')
        myopt.Fourier_mode_when_unstable = true;
    end


    [n,~,N] = size(N0);
    dt = 1/N; % Timestep

    if Type == "White"; Gam = 0; end
    if isscalar(Gam); Gam = Gam*ones(n,1); end

    tCp = (0:(N-1))/N;
    tAp = tCp+0.5*dt;

    % Derivative of covariance function
    dC = (circshift(N0,-1,3) - N0) / dt;

    %% Estimate the dynamics
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
        
        if myopt.Fourier_mode_when_unstable
            % Check if A is stable. If not, then impose Fourier-mode truncation
            A = reshape(AA,n,n,[]);
            [z,h_Check] = CheckStability(A);
            if z == 1 
                myopt.FM = 2;
                AA = SolveEqwithFM(kron(XX',eye(n)),vec(ZZ),myopt.FM);
                h.message = strcat( "Fourier-mode truncation is imposed due to ", h_Check.message );
                h.message_short = h_Check.message_short;
                % disp(h.message)
            end
        end


    end

    A = reshape(AA,n,n,[]);
    
    % Check if the dynamcis is stable.
    if CheckStability_Monodromy(A) == 1

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
            [BT,h_PerSol] = Periodic_Sol( @(t) kron( -Gaminv + At(t), speye(n) ) , @(t) vec( Gaminv ) , tCp );
            BT = reshape(BT,n,n,[]); % This is B^T!!
            B = pagetranspose(BT); % Transpose back to be B!
            h.B = B;

            if h_PerSol.flag == 1 || sum( isnan(B(:)) ) > 0 
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


% Check if the dynamics is stable or not via Floquet exponents
function [z,h] = CheckStability_Monodromy(x)
    
    if max(abs(eig(Monodromy(x)))) > 1 % Unstable dynamics
        z = 1; % Unstable
        h.message = "unstable Floquet multiplier.";
        h.message_short = "US"; % Unstable
    else
        z = 0;
        h = struct();
    end
    
end


% Check if the dynamics is stable entrywise (for the largest n variables) via Frequency analysis
function [z,h] = CheckStability_HighFrequency(x)

    [n,~,N] = size(x);

    x = reshape(x,[],N);
    [~,b] = maxk(vecnorm(x,2,2),n);

    for i = b(:)'
        y = vec(x(i,:));
        y = y - mean(y); % Remove mean

        % High-frequency signal?
        energy = abs(fft(vec(y))).^2;
        totalEnergy = sum(energy);

        onesideEnergy = energy(1:N/2+1);
        onesideEnergy(2:end-1) = 2*onesideEnergy(2:end-1);
        
        EnergyDistribution = onesideEnergy / totalEnergy;
        highFrequencyRatio = sum( EnergyDistribution(end-1:end) );

        if highFrequencyRatio > 0.4
            % EnergyDistribution
            z = 1; % Unstable
            h.message = "excessive high-frequency energy.";
            h.message_short = "HF"; % High frequency
            
            % close all;
            % plot(y)
            % drawnow
            return
        end

    end

    z = 0;
    h = struct();

end

% Check if the dynamics is stable entrywise (for the largest n variables) via peak analysis
function [z,h] = CheckStability_AbnormalPeak(x)

    [n,~,N] = size(x);

    x = reshape(x,[],N);
    [~,b] = maxk(vecnorm(x,2,2),n);

    for i = b(:)'
        y = vec(x(i,:));
        y = y - mean(y); % Remove mean

        % An abnormal peak?
        score = (y).^2;
        totalScore = sum( score );

        ScoreDistribution = score / totalScore;
        PeakRatio1 = sum(maxk(ScoreDistribution,1));
        PeakRatio2 = sum(maxk(ScoreDistribution,2));

        if PeakRatio1 > 0.45 || PeakRatio2 > 0.65
            z = 1; % Unstable
            h.message = "abnormal local curvature.";
            h.message_short = "PK"; % peak
            
            % hold off;
            % yyaxis left;
            % plot(y); hold on;
            % hold off;
            % yyaxis right
            % plot(ScoreDistribution); hold on;
            % ylim([0,0.7])
            % title(strcat(num2str(i)," --- ",num2str(100*PeakRatio2,'%3.1f')))
            % 
            % drawnow
            % pause(0.2)
            
            return
        end

    end

    z = 0;
    h = struct();

end


function [z,h] = CheckStability(x)

    % z: 1/0 (unstable/stable)
    % h.messgage (h.message_short): reason for instability

    % Check if the dynamics is stable or not via Floquet exponents
    [z,h] = CheckStability_Monodromy(x);
    if z == 1; return; end

    % Check if the dynamics is stable entrywise (for the largest n variables)
    [z,h] = CheckStability_HighFrequency(x);
    if z == 1; return; end
    
    % Check if the dynamics is stable entrywise (for the largest n variables)
    [z,h] = CheckStability_AbnormalPeak(x);
    if z == 1; return; end
       
    % if the system is stable
    z = 0;
    h = struct();
    h.message = [];
    h.message_short = [];


end





