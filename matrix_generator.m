n=5;
KER_dim = 1;
%% Creating a regular matrix D with specific eigenvalues > 0
min_eig = 1;
max_eig = 10;
Deig = (max_eig - min_eig)*rand(n,1) + min_eig;
Deig = max(Deig,min_eig);
save('smalleig.mat', 'Deig')
D = diag(Deig);
A = rand(n);
[Q,~] = qr(A);
D = Q'*D*Q;
save('regsmall.mat', 'D')


% %% Creating a singular matrix D with specific eigenvalues >= 0
% min_eig = 1;
% max_eig = 10;
% Seig = (max_eig - min_eig)*rand(n-KER_dim,1) + min_eig;
% Seig = max(Seig,min_eig);
% Seig = [Seig;zeros(KER_dim, 1)];
% save('S3eig.mat', 'Seig')
% S = diag(Seig);
% A = rand(500);
% [Q,~] = qr(A);
% S = Q'*S*Q;
% save('singmatS3.mat', 'S')

