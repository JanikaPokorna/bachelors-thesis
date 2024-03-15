n=500;
KER_dim = 10;
rho = 0.80;
a = 0.1;
b = 100;
%% Creating a regular matrix D with specific eigenvalues > 0
% min_eig = 1;
% max_eig = 10;
% % Deig = (max_eig - min_eig)*rand(n,1) + min_eig;
% % Deig = max(Deig,min_eig);
% D = strakos(n,a,b,rho)
% save('diagonal_matrix_with_eigD5.mat', 'D')
% A = rand(n);
% [Q,~] = qr(A);
% D = Q'*D*Q;
% save('regsmallD5.mat', 'D')

%% Creating a singular matrix D with specific eigenvalues >= 0
min_eig = 1;
max_eig = 10;
% Seig = (max_eig - min_eig)*rand(n-KER_dim,1) + min_eig;
% Seig = max(Seig,min_eig);
% Seig = [Seig;zeros(KER_dim, 1)];
S = strakos(n-KER_dim,a,b,rho);
result_matrix = zeros(n);
result_matrix(1:n-KER_dim, 1:n-KER_dim) = S;
save('diagonal_matrix_with_eig.mat', 'result_matrix')
A = rand(n);
[Q,~] = qr(A);
spanA = Q(:,1:n-KER_dim);
kerA = Q(:,n-KER_dim:n);
S = Q'*result_matrix*Q;
save('singmatS6.mat', 'S')
save('spanS6.mat','spanA')
save('kerS6.mat','kerA')

