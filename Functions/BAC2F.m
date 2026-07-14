function F = BAC2F(B,A,C,g)

    if nargin == 3
        F = @(t,X) f(B(t),X) + A(t)*X + C(t);
    elseif nargin == 4
        F = @(t,X) f(B(t),X) + A(t)*X + C(t) + g(t,X);
    else 
        disp('Check input argument!')
    end

    
end


function y = f(B,x)
    n = length(x);
    y = zeros(n,1);
    for i = 1:n
        for j = 1:n
            for k = 1:n
                y(i) = y(i) + B(i,j,k)*x(j)*x(k);
            end
        end
    end
end