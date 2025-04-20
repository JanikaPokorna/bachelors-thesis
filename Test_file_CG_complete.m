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
tol = 1e-10;

A = gallery('neumann',sizeA);

[V,D] = eigs(A,sizeA);

V = [V(:,end),V(:,1:end-1)];    
[Q,~] = qr(V);      % ortogonalizace proti budoucimu kerA (=posledni vl.vektor)
ker_A = Q(:,1);
span_A = Q(:,2:end);
A = span_A*D(1:end-1,1:end-1)*span_A';

%right hand side vector b
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
new_matrix = mmread('matrices/matrix_market/positive-definite/s2rmt3m1.mtx');
runName = 'Experiment_s2rmt3m1';
maxiter = 20000;
betas = [0, 1e-2];


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

%% K cemu to vubec konverguje? nedava smysl, nejde porovnavat takto ig
n = 500;
x0 = rand(n,1);
ker_dim = 1;
maxiter = 550;

runName = 'Experiment_strakos_nonconvergent';

rho = 0.8; %the smaller this is, the more eigenval are close, should be below 1
a = 5;
c = 100;
[A,D,spanA,kerA] = singular_strakos(n,ker_dim,a,c,rho); % creates strakos matrix with ker of dimension 1


%create right-hand side vector b
b = make_multi_vector_b(spanA,kerA,betas,1);
figure;
semilogy(diag(D), 'or');
grid on;

run_CG_nonconvergent(runName,A,b,x0,maxiter,tol,betas,spanA,kerA)
%% strakos s vetsi dimenzi jadra
n = 500;
x0 = zeros(n,1);
ker_dim = 3;
maxiter = 550;
tol = 1e-8;

runName = 'Experiment_strakos_larger_kernel';

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
%% Orthodir
betas = [0, 1e-4, 1e-3];
n = 500;
x0 = zeros(n,1);
ker_dim = 10;
maxiter = 550;

runName = 'Experiment_Orthodir';

rho = 0.8;
a = 5;
c = 100;
[A,D,spanA,kerA] = singular_strakos(n,ker_dim,a,c,rho);


%create right-hand side vector b
b = make_multi_vector_b(spanA,kerA,betas);


run_analysis_Orthodir(runName,A,b,x0,maxiter,tol,betas,spanA,kerA)

%% Orthodir pro Neumanna
runName = 'Experiment_neumann_orthodir';
n = 25;
sizeA = n^2;
x0 = zeros(sizeA,1);
maxiter = 500;
tol = 1e-10;

A = gallery('neumann',sizeA);

[V,D] = eigs(A,sizeA);

V = [V(:,end),V(:,1:end-1)];    
[Q,~] = qr(V);      
% ortogonalizace proti budoucimu kerA (=posledni vl.vektor)
ker_A = Q(:,1);
span_A = Q(:,2:end);
A = span_A*D(1:end-1,1:end-1)*span_A';

b = make_multi_vector_b(span_A',ker_A',betas);

run_analysis_Orthodir(runName,A,b,x0,maxiter,tol,betas,span_A',ker_A')
%% Orthodir for non-singular matrices
n = 500;
x0 = zeros(n,1);
ker_dim = 0;
maxiter = 1000;



rho = 0.8; %the smaller this is, the more eigenval are close, should be below 1
a = 5;
c = 100;
[A,D,spanA,kerA] = singular_strakos(n,ker_dim,a,c,rho); % creates strakos matrix with ker of dimension 1


%create right-hand side vector b
b = make_multi_vector_b(spanA,kerA,betas);
figure;
semilogy(diag(D), 'or');
grid on;
% runName = 'Experiment_strakos_nonsing_CG';
% run_analysis_CG(runName,A,b,x0,maxiter,tol,betas,spanA,kerA)
runName = 'Experiment_strakos_nonsing_orthodir';
run_analysis_Orthodir(runName,A,b,x0,maxiter,tol,betas,spanA,kerA)
