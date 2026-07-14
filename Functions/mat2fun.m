function Xt = mat2fun(X,t)

    % Make the matrix X a function along the 3rd dimension: X(i,j,k) -> X(i,j,t)

    t2idx = @(t) min([floor(mymod(100*t+1,100)),100]); % "mod" may suffer from numerical floating point error. e.g., "mymod(100*-5.5511e-17+1,100) = 101".
    D_fine = myinterp1( t, X, 0.01*(0:99) );
    Xt = @(t) D_fine(:,:,t2idx(t));

end