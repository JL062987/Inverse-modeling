function [K,h] = CS_CorrFcn(Type,dt,maxlag,A,Q,D,Gam,tAp)

    if isempty(Gam)
        Gam = 0;
    end

    if isa(A,'function_handle')
        [n,~] = size(A(1));
    else
        [n,~] = size(A);
    end
    N = length(tAp);

    idx = 1:n;

    if Type == "White" || Gam == 0
        C = Cov(A,D,tAp);
    elseif Type == "Colored" || Type == "CW"
        if isempty(D) || Type == "Colored"; D = []; end 
        [A,D] = CS_to_Aug_System(Gam,A,Q,D);
        C = Cov(A,D,tAp);
    end


    if maxlag == 0
        G = 0;
        K = C;
    else
        G = Green(A,tAp,dt,maxlag);
        K = zeros(size(G));
        for i = 1:N
            K(:,:,i,:) = pagemtimes(G(:,:,i,:),C(:,:,i));
        end
    end

    K = K(idx,idx,:,:);

    h = struct();
    h.C = C;
    h.G = G;
       

end

function C = Cov(A,Q,tAp)

    [n,~] = size(A(1));
    N = length(tAp);
    I = speye(n);

    L = @(t) kron(I,A(t)) + kron(A(t),I);
    c = @(t) +2*vec(Q(t));

    % tic
    [~,xp] = ode45(@(t,y) L(t)*y + c(t), [-1,25+tAp], vec(eye(n)));
    % disp(num2str(toc))
    xp = xp';
    xp = xp(:,2:end);
    C = reshape(xp,n,n,N);
    

end

function G = Green(A,tAp,dt,maxlag)
    
    [n,~] = size(A(0));
    N = length(tAp);
    L = zeros(n,n,N);

    for i = 1:N
        t0 = tAp(i); % Initial time
        M = @(t) kron(eye(n),A(t));
        [~,xp] = ode45(@(t,y) M(t)*y, [t0,t0+dt], vec(eye(n)));
        xp = xp';
        L(:,:,i) = reshape(xp(:,end),n,n);
    end


    G = zeros(n,n,N,maxlag+1);
    G(:,:,:,1) = repmat(eye(n),1,1,N);

    for i = 1:N
        for j = 1:maxlag
            G(:,:,i,j+1) = L(:,:,mymod(i+j-1,N))*G(:,:,i,j);     
        end
    end

end