function [G,Gd,Gs] = CS_Evolution_Operator(Type,C,Coupling,tC,A,tA,Q,D,Gam,tlead,tphase)


    % tlead must be evenly distributed!
    
    [n,~,N] = size(A);
    
    if isscalar(Gam)
        Gam = Gam*ones(n,1);
    end
    
    Gd = zeros(n,n,length(tlead),length(tphase));
    Gs = zeros(n,n,length(tlead),length(tphase));
    
    if Type == "White"
        for i = 1:length(tphase) % For each theta
        
            M = 8;
        
            t = tphase(i) + linspace(0,tlead(end),M*length(tlead)+1); % Make the integration timestep smaller
            dt = t(2) - t(1);
        
            Aq = myinterp1(tA,A,mod(t,1));
            
            EvoMatrix = eye(n);
            for k = 1:(M*length(tlead))
                EvoMatrix = expm(Aq(:,:,k)*dt) * EvoMatrix;
                if mod(k,M) == 0
                    j = myquotient(k,M);
                    Gd(:,:,j,i) = EvoMatrix;
                end
            end
        end
    
    elseif Type == "Colored" || Type == "CW" 
    
        L = zeros(2*n,2*n,N);
    
        sq2Q = sqrtm(2*Q);
    
        idx1 = 1:n;
        idx2 = (n+1):(2*n);
        
        L(idx1,idx1,:) = A;
        L(idx1,idx2,:) = repmat(sq2Q,1,1,N);
        L(idx2,idx2,:) = repmat(diag(-1./Gam),1,1,N);
    
        for i = 1:length(tphase)
        
            M = 8;
        
            t = tphase(i) + linspace(0,tlead(end),M*length(tlead)+1);
            dt = t(2) - t(1);
        
            Lq = myinterp1(tA,L,mod(t,1));
            Cq = myinterp1(tC,C,mod(t,1));
            Pq = myinterp1(tC,Coupling,mod(t,1));
            
            EvoMatrix = eye(2*n);
            for k = 1:(M*length(tlead))
                EvoMatrix = expm(Lq(:,:,k)*dt) * EvoMatrix;
                if mod(k,M) == 0
                    j = myquotient(k,M);
                    Gd(:,:,j,i) = EvoMatrix(idx1,idx1);
                    Gs(:,:,j,i) = EvoMatrix(idx1,idx2)*(Pq(:,:,1)/Cq(:,:,1));
                end
            end
        end
        
    end
    
    G = Gd + Gs;

end