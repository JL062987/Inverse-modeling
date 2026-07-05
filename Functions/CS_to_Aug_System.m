function [L,S,sq2S] = CS_to_Aug_System(Gam,A,Q,D)

    % Input and output have the same format: matrix -> matrix, function -> function

    if isa(A,"function_handle")

        [n,~,~] = size(A(1));

        if nargin < 4 || isempty(D) 
            D = @(t) zeros(n); sq2D = @(t) zeros(n);
        else
            sq2D = @(t) sqrtm(2*D(t));
        end

        sq2Q = sqrtm(2*Q); Gaminv = (1/Gam)*eye(n);
        L = @(t) [A(t) sq2Q; zeros(n), -Gaminv ];

        Gaminvsq = 1/Gam.^2;
        S = @(t) [ D(t) zeros(n); zeros(n), 0.5*Gaminvsq*eye(n) ];

        sq2S = @(t) [ sq2D(t) zeros(n); zeros(n) Gaminv ];


    else

        [n,~,N] = size(A);
        idx1 = 1:n; idx2 = (n+1):(2*n);
    
    
        if nargin < 4; D = zeros(n,n,N); end
    
        L = zeros(2*n,2*n,N);
        L(idx1,idx1,:) = A;
        L(idx1,idx2,:) = repmat( sqrtm(2*Q), 1, 1, N );
        L(idx2,idx2,:) = repmat( (-1./Gam)*eye(n) , 1, 1, N );
        
        S = zeros(2*n,2*n,N);
        S(idx1,idx1,:) = D;
        S(idx2,idx2,:) = repmat((1./(2*Gam.^2))*eye(n),1,1,N);
    
        sq2S = [];
        % if norm(D(:)) < 1e-10
        %     sq2S = [ zeros(n), zeros(n);   zeros(n), diag( 1./Gam )];
        % else
        %     sq2S = [ sqrtm(2*D), zeros(n);   zeros(n), diag( 1./Gam )];
        % end

    end

end