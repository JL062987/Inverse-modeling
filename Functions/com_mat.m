function P = com_mat(m, n)

% determine permutation applied by K
A = reshape(1:m*n, m, n);
v = reshape(A', 1, []);

% apply this permutation to the rows (i.e. to each column) of identity matrix
P = eye(m*n);
P = P(v,:);

end