%% test_file strakos singular multiple plots

n = 500;
ker_dim = 1;
maxiter = 550;
tol = 1e-15;
x0 = zeros(n,1);
rho = 0.8; %the smaller this is, the more eigenval are close, should be below 1
a = 5;
c = 100;
deltas = [0];
[S,D,spanA,kerA] = singular_strakos(n,ker_dim,a,c,rho); % creates strakos matrix with ker of dimension 1
diagonal_values = diag(D);

% create right-hand side vector b
b = make_multi_vector_b(spanA,kerA,deltas);


% Plot the eigen values of symmetric PD strakos matrix 
figure(2);
semilogy(diag(D), 'or');
title('Diagonal values of the matrix');
xlabel('Index');
ylabel('Diagonal Value');
grid on;


% Do CG for each right-hand side vector b(i)
X_matrices = cell(size(b,2), 1);
Gamma_matrices = cell(size(b,2), 1);
for i = 1:size(b,2)
    [x, X_i,l,~,~,Gamma_i] = conjugate_grad(S, b(:,i),x0,maxiter,tol);
    if i == 1
        converged_x = x;
        len = l;
    end
    X_matrices{i} = X_i;
    Gamma_matrices{i,1} = Gamma_i;
end

% Make error matrices
Error_matrices = cell(size(b,2), 1);
for k = 1:size(b,2)
    error_matrix = zeros(1,maxiter);
    X = X_matrices{k};
    for i = 1:maxiter
        A_norm_xi = sqrt((converged_x - X(:,i))'*S*(converged_x - X(:,i)));
        error_matrix(1,i) = A_norm_xi/(sqrt((converged_x - x0)'*S*(converged_x - x0)));
        if (i > 5 && error_matrix(1,i)==1)
            error_matrix(1,i:end) =0;
            break
        end
    end
    Error_matrices{k} = error_matrix;
end

colors = {'r', 'g', 'b', 'm', 'c'};

%check dimensions
size(Error_matrices{1}(1:maxiter))
size(1:maxiter)

figure(3)
hold on;
for i = 1:size(b,2)
    semilogy(1:maxiter, Error_matrices{i}(1:maxiter),['o-', colors{i}]);
    hold on;
end
title('Relative error for different delta values');
xlabel('Step k');
ylabel('||x - x_i||_A / ||x - x_0||_A');
ylim([0, 8]);
xlim([0,len])
legend('\delta = 0', '\delta = 1e-2', '\delta = 1e-4','\delta = 1e-6', 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;

figure(4)
hold on;
for i = 1:size(b,2)
    plot(1:len+1, Gamma_matrices{i}(1:len+1),['o-', colors{i}]);
    hold on;
end
xlabel('Step k')
ylabel('Gamma values for different delta values')
title('Values of Gamma')
legend('\delta = 0', '\delta = 1e-2', '\delta = 1e-4','\delta = 1e-6', 'Location', 'best');
hold off;
