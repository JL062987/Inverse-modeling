function [y, lambda] = yeo_johnson_transform(x)

    x = x(:);
    N = length(x);
    jacobian_sum = sum( sign(x) .* log(abs(x) + 1) );

    % Check if the input contains all zeros (cannot be transformed)
    if var(x) == 0
        error('Data cannot be constant (variance must be greater than 0).');
    end
    
    % Initialize optimal lambda using maximum likelihood estimation (MLE)
    lambda_range = -2:0.01:2; % Range of lambda to test
    logLikelihood = zeros(size(lambda_range));

    
    for i = 1:length(lambda_range)

        lambda = lambda_range(i);
        y0 = apply_yeo_johnson(x, lambda); % Apply the Yeo-Johnson transformation

        % 1. Calculate MLE variance (using N as denominator instead of N-1)
        var_mle = var(y0,0,1);
        
        % 2. Calculate the log of the Jacobian determinant
        log_jacobian = (lambda - 1) * jacobian_sum;

        % 3. Combine them into the correct profile log-likelihood
        logLikelihood(i) = -0.5 * N * log(var_mle) + log_jacobian;

    end

    % Find the optimal lambda corresponding to the maximum log-likelihood
    [~, idx] = max(logLikelihood);
    lambda = lambda_range(idx);

    % Apply the Yeo-Johnson transformation with the optimal lambda
    y = apply_yeo_johnson(x, lambda);

end


function y = apply_yeo_johnson(x, lambda)

    % Completely vectorized Yeo-Johnson transformation
    y = zeros(size(x));
    
    % Logical masks for positive and negative values
    pos = (x >= 0);
    neg = ~pos;
    
    % Transformation for x >= 0
    if lambda ~= 0
        y(pos) = ((x(pos) + 1).^lambda - 1) / lambda;
    else
        y(pos) = log(x(pos) + 1);
    end
    
    % Transformation for x < 0
    if lambda ~= 2
        y(neg) = -(((-x(neg) + 1).^(2 - lambda) - 1) / (2 - lambda));
    else
        y(neg) = -log(-x(neg) + 1);
    end

end