function [L,S,sq2S] = ST_to_Aug_System(Gam,A,Q,D)

    [n,~] = size(A);

    if isscalar(Gam); Gam = Gam*ones(n,1); end
    if nargin < 4; D = zeros(n); end

    L = [ A, sqrtm(2*Q); zeros(n), -diag(1./Gam) ]; % Augmented dynamical matrix
    S = [ D, zeros(n);   zeros(n), diag( 1./(2*Gam.^2) )]; % Augmented stochastic matrix
    
    if norm(D(:)) < 1e-10
        sq2S = [ zeros(n), zeros(n);   zeros(n), diag( 1./Gam )];
    else
        sq2S = [ sqrtm(2*D), zeros(n);   zeros(n), diag( 1./Gam )];
    end
end