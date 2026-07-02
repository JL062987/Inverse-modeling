function K = CorrFcn_ST_Obs(xvec,maxlag)

    % xvec: var * time (e.g. monthly data)
    % The lag variable is sampled as the same as input xvec.
    
    [n,~] = size(xvec);
    K = zeros(n,n,2*maxlag+1); % Observed correlation
    
    for i = 1:n
        for j = 1:n
            xi = xvec(i,:);
            xj = xvec(j,:);
            %xi = zscore(xi);
            %xj = zscore(xj);
            K(i,j,:) = xcorr(xi,xj,maxlag,"unbiased");
        end
    end

    K = K(:,:,maxlag+1:end);

end