function [ Atilde, btilde, x2_comp ] = Schur_Complement(A,b,k)

    [n1,n2] = size(A);

    if n1 == n2
        n = n1;
    else
        disp('Dimension error: Schur complement.')
    end
    
    i1 = 1:k;
    i2 = (k+1):n;

    A11 = A(i1,i1);
    A12 = A(i1,i2);
    A21 = A(i2,i1);
    A22 = A(i2,i2);
    b1 = b(i1); 
    b2 = b(i2);

    Atilde = A11 - A12 * ( A22 \ A21 ) ;
    btilde = b1 - A12 * ( A22 \ b2 ) ;

    x2_comp = @(x1) A22 \ ( b2 - A21 * x1 );

end