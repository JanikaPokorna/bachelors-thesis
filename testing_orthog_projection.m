%% projecting vector onto subspaces using pseudoinverse
n = 100;
m = 80;
x = rand(n,1);
A = eye(n);
spanA = A(:,1:m);
kerA = A(:,m+1:n);

P = spanA * pinv(spanA);
Q = kerA * pinv(kerA);
% Project x onto the subspace spanned by the basis vectors
spanA_x = P * x;
kerA_x = Q * x;

reconstructed_x = spanA_x + kerA_x;

% Check if the original x equals the reconstructed x
disp(all(x == reconstructed_x));