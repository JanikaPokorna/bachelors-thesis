%% test_file strakos singular multiple plots

n = 500;
ker_dim = 1;
x0 = zeros(n,1);
rho = 0.8; %the smaller this is, the more eigenval are close, should be below 1
a = 5;
c = 100;
deltas = [0,1e-2,1e-4,1e-6];
[S,D,spanA,kerA] = singular_strakos(n,ker_dim,a,c,rho); % creates strakos matrix with ker of dimension 1
diagonal_values = diag(D);

% create right-hand side vector b
b = zeros(n,size(deltas,2));
for i = 1:size(deltas,2)
    b(:,i) = make_vector_b(spanA,kerA,deltas(i));
end

eigen_values = eig(S);

% Plot the eigen values of symmetric PD strakos matrix 
figure(2);
semilogy(eigen_values, 'or');
title('Diagonal values of the matrix');
xlabel('Index');
ylabel('Diagonal Value');
grid on;


% Do CG for each right-hand side vector b(i)
X_matrices = cell(size(b,2), 1);
Gamma_matrices = cell(size(b,2), 1);
for i = 1:size(b,2)
    [x, X_i,l,~,Gamma_i] = conjugate_grad(S, b(:,i));
    if i == 1
        converged_x = x;
        len = l;
    end
    X_matrices{i} = X_i;
    Gamma_matrices{i} = Gamma_i;
end

% Make error matrices
Error_matrices = cell(size(b,2), 1);
for k = 1:size(b,2)
    error_matrix = zeros(1,len);
    X = X_matrices{k};
    for i = 1:len
        A_norm_xi = sqrt((converged_x - X(:,i))'*S*(converged_x - X(:,i)));
        error_matrix(1,i) = A_norm_xi/(sqrt((converged_x - x0)'*S*(converged_x - x0)));
    end
    Error_matrices{k} = error_matrix;
end

colors = {'r', 'g', 'b', 'm', 'c'};

figure(3)
hold on;
for i = 1:size(b,2)
    semilogy(1:len, Error_matrices{i}(1:len),['o-', colors{i}]);
    hold on;
end
title('Relatove error for different delta values');
xlabel('Step k');
ylabel('||x - x_i||_A / ||x - x_0||_A');
ylim([0, 5]);
xlim([0,50])
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
hold off;
