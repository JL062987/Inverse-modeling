function K = Phase_Averaging(K,myopt)

    [n,~,N,M] = size(K);

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


    % Smooth the correlation function
    lowpassFilter = designfilt('lowpassiir', 'FilterOrder', 8, ...
        'HalfPowerFrequency', fpass, 'SampleRate', fs);
    
    % Lowpass filter on the observed corelation function in the phase direction
    for k = 1:M
        for i = 1:n
            for j = 1:n
                x = vec(K(i,j,:,k));
    
                y = filtfilt(lowpassFilter, [x;x;x]);
                y = vec(y(N+1:(2*N)));
        
                K(i,j,:,k) = y;
            end
        end
    end

end