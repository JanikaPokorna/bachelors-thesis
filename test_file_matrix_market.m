%% test file matrix market - bcsstk04

betas = [0, 1e-6, 1e-4, 20];
tol = 1e-15;
maxiter = 3000;


%% load matrix

% new_data = importdata('matrices/matrix_market/positive-definite/bcsstk04.mat');
% new_matrix = new_data.A;
% n = size(new_matrix,1);

new_matrix = mmread('matrices/matrix_market/positive-definite/bcsstm26.mtx');
% new_matrix = mmread('matrices/matrix_market/positive-definite/1138_bus.mtx');
% new_matrix = mmread('matrices/matrix_market/positive-definite/nos7.mtx');
% new_matrix = mmread('matrices/matrix_market/positive-definite/s2rmt3m1.mtx');

n = size(new_matrix,1);
x0 = zeros(n,1);


%% find eigenvalues

[V,D] = eigs(new_matrix,n);

%% modify to semi-definite
D(n,n) = 0;
M = V * D * V';

%% create right hand side

span_M = V(:,1:n-1); % sloupce
ker_M = V(:,n);
% is_in_kernel = norm(M * ker_M) < 1e-10

b = make_multi_vector_b(span_M',ker_M',betas); % musim do funkce dosadit radky
b_ker = ker_M*ker_M' * b; %VSUDE NULY, ALE ASI FUNGUJE, JEN JSOU ROZDILY MINIMALNI
%% plot of eigenvalues

figure(1);
semilogy(diag(D), 'or');
title('Diagonal values of the matrix');
xlabel('Index');
ylabel('Diagonal Value');
grid on;

%% CG
X_matrices = cell(size(b,2), 1);
Gamma_matrices = cell(size(b,2), 1);
P_matrices = cell(size(b,2), 1);
P_projected_ker = cell(size(b,2), 1);
X_projected_ker = cell(size(b,2), 1);
R_matrices = cell(size(b,2), 1);
R_projected_span = cell(size(b,2), 1);


for i = 1:size(b,2)
    [x, X_i,l,P_i,R_i,Gamma_i] = conjugate_grad(M, b(:,i),x0,maxiter,tol);
    if i == 1
        converged_x = x;
        len = l;
    end
    X_matrices{i} = X_i;
    Gamma_matrices{i,1} = Gamma_i;
    P_matrices{i,1} = P_i;
    R_matrices{i,1} = R_i;
    % Oddělení komponenty v jádře
    column_norms = norm(ker_M', 2);  % normy sloupců
    is_orthonormal_cols = all(abs(column_norms - 1) < eps);
    if is_orthonormal_cols == true
        %kontrola dimenzí
        size(ker_M')
        size(X_matrices{i}) 
        X_projected_ker{i} = (ker_M*ker_M')* X_matrices{i};
        P_projected_ker{i} = (ker_M*ker_M')*P_matrices{i};
        R_projected_span{i} = (span_M*span_M')*R_matrices{i};
    end
end

%% error matrices

Error_matrices = cell(size(b,2), 1);
for k = 1:size(b,2)
    error_matrix = zeros(1,maxiter);
    X = X_matrices{k};
    for i = 1:maxiter
        M_norm_xi = sqrt((converged_x - X(:,i))'*M*(converged_x - X(:,i)));
        error_matrix(1,i) = M_norm_xi/(sqrt((converged_x - x0)'*M*(converged_x - x0)));
        if (i > 5 && error_matrix(1,i)==1)
            error_matrix(1,i:end) = 0;
            break
        end
    end
    Error_matrices{k} = error_matrix;
end

%Error matrices for kernel component of x_k
Error_matrices_kernel = cell(size(b,2), 1);
for k = 1:size(b,2)
    error_matrix_ker = zeros(1,maxiter);
    X = X_projected_ker{k};
    first_vector_x1 = norm(X(:,2),2);
    for i = 1:maxiter
        Norm_xi = norm(X(:,i),2);
        error_matrix_ker(1,i) = Norm_xi/norm(first_vector_x1,2);
        if (i > 5 && error_matrix_ker(1,i)==1)
            error_matrix_ker(1,i:end) = 0;
            break
        end
    end
    Error_matrices_kernel{k} = error_matrix_ker;
end

%A-Norm of vectors p_k component in kernel (A=S)
Norm_matrices_p = cell(size(b,2), 1);
for k=1:size(b,2)
    P = P_projected_ker{k};
    p_norm_matrix = zeros(1,maxiter);
    size(p_norm_matrix)
    for i = 1:maxiter
        Norm_p = norm(P(:,1),2);
        p_norm_matrix(1,i) = Norm_p;
    end
Norm_matrices_p{k} = p_norm_matrix;
end

%comparing Gamma components
Span_components = cell(size(b,2), 1);
Ker_components = cell(size(b,2), 1);
for k = 1: size(b,2)
    span_component_matrix = zeros(1,maxiter);
    ker_component_matrix = zeros(1,maxiter);
    R = R_matrices{k};
    P = P_matrices{k};
    b_kernel = b_ker(k);
    for i=1:maxiter
        denominator = P(:,i)'*M*P(:,i);
        span_component_matrix(1,i) = (R(:,i)'*R(:,i))/denominator;
        ker_component_matrix(1,i) = (b_kernel*b_kernel')/denominator;
    end
    Span_components{k} = span_component_matrix;
    Ker_components{k} = ker_component_matrix;
end

colors = {'r', 'g', 'b', 'm', 'c'};

%check dimensions
size(Error_matrices{1}(1:maxiter));
size(1:maxiter);

%% plotting 

figure(2)
hold on;
for i = 1:size(b,2)
    semilogy(1:maxiter, Error_matrices{i}(1:maxiter),['o-', colors{i}]);
    hold on;
end
title('Relative error for different beta values');
xlabel('Step k');
ylabel('||x - x_i||_A / ||x - x_0||_A');
ylim([0, inf]);
xlim([0,len]);
legend('\beta = 0', '\beta = 1e-6', '\beta = 1e-4','\beta = 1e-2', 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;

figure(3)
hold on;
for i = 1:size(b,2)
    plot(1:len+1, Gamma_matrices{i}(1:len+1),['o-', colors{i}]);
    hold on;
end
xlabel('Step k')
ylabel('Gamma values')
title('Values of Gamma')
legend('\beta = 0', '\beta = 1e-6', '\beta = 1e-4','\beta = 1e-2', 'Location', 'best');
hold off;

figure(4)
hold on;
for i = 1:size(b,2)
    semilogy(1:maxiter, Error_matrices_kernel{i}(1:maxiter),['o-', colors{i}]);
    hold on;
end
title('Relative error - kernel components of x_i');
xlabel('Step k');
ylabel('||x_i||/ ||x_0||');
ylim([0, inf]);
xlim([0,len+30])
legend('\beta = 0', '\beta = 1e-6', '\beta = 1e-4','\beta = 1e-2', 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;

for i = 1:size(b,2)
    figure(6+i);
    hold on;
    h_kernel = semilogy(1:maxiter, Ker_components{i}(1:maxiter),'.-r');
    h_span = semilogy(1:maxiter, Span_components{i}(1:maxiter),'.-b');
    title('Comparison of \gamma components');
    xlabel('Step k', 'FontName', 'AvantGarde');
    ylabel('Value', 'FontName', 'AvantGarde');
    ylim([0, inf]);  
    xlim([0, len]);
    set(gca, 'YScale', 'log');
    % set([ h_kernel h_span], 'LineWidth', 1.5, 'MarkerSize', 6);
    set(h_kernel, 'LineWidth', 1.5, 'MarkerSize', 6);
    set( h_kernel,'MarkerFaceColor',[.0 0 .8],'MarkerEdgeColor','r');
    % set( h_span, 'MarkerFaceColor',[.0 0 .8],'MarkerEdgeColor','b');
    legend('Kernel Component', 'Span Component', 'Location', 'best');
    hold off;
end