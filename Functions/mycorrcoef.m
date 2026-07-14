function z = mycorrcoef(u,v) 

    % Compute correlation coefficient row wise

    Take21 = @(v) v(2,1);

    [n,m] = size(u);

    if n == 1 || m == 1 % a vector
        z = Take21( corrcoef( u(:), v(:), 'Rows', 'Complete' ) );
    
    else    
        z = zeros(n,1);    
        for i = 1:n
            z(i) = Take21( corrcoef( vec(u(i,:)), vec(v(i,:)), 'Rows', 'Complete' ) );
        end
    end

end