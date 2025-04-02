%%Test file shortened script
tol = 1e-15;
betas = [0, 1e-4, 1e-3, 1e-2];
%%
n = 500;
x0 = zeros(n,1);
ker_dim = 1;
maxiter = 550;

runName = 'Experiment_strakos';

rho = 0.8; %the smaller this is, the more eigenval are close, should be below 1
a = 5;
c = 100;
[A,D,spanA,kerA] = singular_strakos(n,ker_dim,a,c,rho); % creates strakos matrix with ker of dimension 1
diagonal_values = diag(D);

%create right-hand side vector b
b = make_multi_vector_b(spanA,kerA,betas);
figure;
semilogy(diag(D), 'or');
grid on;

run_analysis_CG(runName,A,b,x0,maxiter,tol,betas,spanA,kerA)
%%
runName = 'Experiment_neumann';
n = 25;
sizeA = n^2;
x0 = zeros(sizeA,1);
maxiter = 300;

A = gallery('neumann',sizeA);

%right hand side vector b
[V,D] = eigs(A,sizeA);
span_A = V(:,1:sizeA-1);
ker_A = V(:,sizeA);
b = make_multi_vector_b(span_A',ker_A',betas);

figure;
semilogy(diag(D), 'or');
grid on;

run_analysis_CG(runName,A,b,x0,maxiter,tol,betas,span_A',ker_A')
%% 
maxiter = 4000;
% new_matrix = mmread('matrices/matrix_market/positive-definite/bcsstm26.mtx');
% runName = 'Experiment_bcsstm26';
% new_matrix = mmread('matrices/matrix_market/positive-definite/1138_bus.mtx');
% runName = 'Experiment_1138_bus';
% new_matrix = mmread('matrices/matrix_market/positive-definite/nos7.mtx');
% runName = 'Experiment_nos7';
% new_matrix = mmread('matrices/matrix_market/positive-semidefinite/nos4.mtx');
% runName = 'Experiment_nos4';
% maxiter = 300;
% new_matrix = mmread('matrices/matrix_market/positive-definite/s2rmt3m1.mtx');
% runName = 'Experiment_s2rmt3m1';
% maxiter = 20000;
% betas = [0, 1e-2];


n = size(new_matrix,1);
x0 = zeros(n,1);

% modify to semi-definite
[V,D] = eigs(new_matrix,n);
D(n,n) = 0;
A = V * D * V';
% create right hand side
span_A = V(:,1:n-1);
ker_A = V(:,n);
b = make_multi_vector_b(span_A',ker_A',betas);

figure;
semilogy(diag(D), 'or');
grid on;

run_analysis_CG(runName,A,b,x0,maxiter,tol,betas,span_A',ker_A')
%% Semidefinite matrix  / NENI SPSP!!!!!


A = mmread('matrices/matrix_market/positive-semidefinite/nos4.mtx');
runName = 'Experiment_nos4';
sizeA = size(A,1);
n = sizeA;
x0 = zeros(n,1);
maxiter = 300;

%right hand side vector b
[V,D] = eigs(A,sizeA);
span_A = V(:,1:sizeA-1);
ker_A = V(:,sizeA);
b = make_multi_vector_b(span_A',ker_A',betas);

figure;
semilogy(diag(D), 'or');
grid on;

run_analysis_CG(runName,A,b,x0,maxiter,tol,betas,span_A',ker_A')


