function [x_mat,t_vec] = SDE_integrate_color(f, g, tf, x0, h, eta)

    % Integrate x' = f(t, x) + sqrt(2g)*eta using classical RK4 or Euler forward
    %
    %   f     : function handle f(t, x)
    %   tspan : tf or [ti, tf]
    %   x0    : initial condition (column vector)
    %   h     : time step
    %
    %   x_mat : solution, each column corresponds to time in t_vec
    %   t_vec : time grid (row vector)

    if isscalar(tf)
        t_vec = 0:h:tf;
    else
        t_vec = tf(1):h:tf(2);
    end
    
    n_steps = numel(t_vec);
    n_dim   = numel(x0);
    
    x_mat = zeros(n_dim,n_steps);
    x_mat(:,1) = x0(:);

    sq2g = sqrtm(2*g);
    
    for n = 1:n_steps-1
        tn = t_vec(n);
        xn = x_mat(:,n);
        etan = eta(:,n);
    
        % Euler forward
        x_mat(:,n+1) = (xn + h * f(tn, xn) + h*sq2g*etan ).';

    end

end