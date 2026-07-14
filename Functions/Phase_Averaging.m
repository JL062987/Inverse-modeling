function Kout = Phase_Averaging(K,myopt)


    if nargin < 2
        myopt = struct();
    end

    % Setup parameters
    if isfield(myopt,'fs')
        fs = myopt.fs;
    else
        fs = 12; % Assuming monthly data
    end

    if isfield(myopt,'fpass')
        fpass = myopt.fpass;
    else
        fpass = 0.2*fs;
    end


    % Define a filter
    lowpassFilter = designfilt('lowpassiir', 'FilterOrder', 8, ...
        'HalfPowerFrequency', fpass, 'SampleRate', fs);
    
   
    % [n,~,N,M] = size(K);
    % 
    % % Lowpass filter on the observed corelation function in the phase direction
    % for k = 1:M
    %     for i = 1:n
    %         for j = 1:n
    %             x = vec(K(i,j,:,k));
    % 
    %             y = filtfilt(lowpassFilter, [x;x;x]);
    %             y = vec(y(N+1:(2*N)));
    % 
    %             K(i,j,:,k) = y;
    %         end
    %     end
    % end

    % Smooth the correlation function

    % By default, the input matrix is of ivar x ivar x phase x maxlead and phase
    % average is applied along the phase varialble. However, it can be
    % specified by 'dim' via myopt. 


    % Target dimension
    if isfield(myopt,'dim')
        dim = myopt.dim;
    else
        dim = 3; 
    end
    
    nd = ndims(K);

    % Move the phase dimension to the first dimension
    order = [dim, 1:dim-1, dim+1:nd];
    Kperm = permute(K, order);

    % Store the permuted size
    sz = size(Kperm);
    N  = sz(1);

    % Each column is one phase-dependent series
    Kmat = reshape(Kperm, N, []);

    % Periodic extension
    Kmat = repmat(Kmat, 3, 1);

    % Filter all columns along the phase direction
    Kmat = filtfilt(lowpassFilter, Kmat);

    % Retain the middle period
    Kmat = Kmat(N+1:2*N, :);

    % Restore the permuted array shape
    Kperm = reshape(Kmat, sz);

    % Restore the original dimension order
    Kout = ipermute(Kperm, order);

end