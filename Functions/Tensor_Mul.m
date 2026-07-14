function [M,x] = Tensor_Mul(A,B,k)

    % Convert A x_k B into [L_B^k]*vec(A).

    nA = ndims(A);
    nB = ndims(B);
    szA = size(A);
    szB = size(B);
    
    M = kron( ( reshape(B,[prod(szB(1:k)),prod(szB(k+1:nB))]) )', eye(prod(szA(1:(nA-k)))) );
    x = vec(A);

end