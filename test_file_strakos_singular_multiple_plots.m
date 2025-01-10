%% test_file strakos singular multiple plots

n = 500;
ker_dim = 1;
maxiter = 550;
tol = 1e-15;
x0 = zeros(n,1);
rho = 0.8; %the smaller this is, the more eigenval are close, should be below 1
a = 5;
c = 100;
betas = [0, 1e-6, 1e-4, 1e-2];
[S,D,spanA,kerA] = singular_strakos(n,ker_dim,a,c,rho); % creates strakos matrix with ker of dimension 1
diagonal_values = diag(D);

%% create right-hand side vector b
b = make_multi_vector_b(spanA,kerA,betas);
% b kernel component
b_ker = kerA'*kerA * b;
% b_kernel'*b_kernel;

%% Plot the eigenvalues of symmetric PD strakos matrix 
figure(2);
semilogy(diag(D), 'or');
title('Diagonal values of the matrix');
xlabel('Index');
ylabel('Diagonal Value');
grid on;


%% Do CG for each right-hand side vector b(i)
% Zavedení proměnných
X_matrices = cell(size(b,2), 1);
Gamma_matrices = cell(size(b,2), 1);
Delta_matrices = cell(size(b,2), 1);
P_matrices = cell(size(b,2), 1);
P_projected_ker = cell(size(b,2), 1);
R_projected_span = cell(size(b,2), 1);
X_projected_ker = cell(size(b,2), 1);
X_projected_span = cell(size(b,2), 1);
R_matrices = cell(size(b,2), 1);
%CG + zaznamenání do proměnných
for i = 1:size(b,2)
    [x, X_i,l,P_i,R_i,Gamma_i,Delta_i] = conjugate_grad(S, b(:,i),x0,maxiter,tol);
    if i == 1
        converged_x = x;
        len = l;
    end
    X_matrices{i} = X_i;
    Gamma_matrices{i,1} = Gamma_i;
    Delta_matrices{i,1} = Delta_i;
    P_matrices{i,1} = P_i;
    R_matrices{i,1} = R_i;
    % Oddělení komponnty v jádře
    column_norms = norm(kerA, 2);  % normy sloupců
    is_orthonormal_cols = all(abs(column_norms - 1) < eps);
    if is_orthonormal_cols == true
        %kontrola dimenzí
        size(kerA)
        size(X_matrices{i}) 
        % X_projected_ker{i} = kerA* X_matrices{i};
        % P_projected_ker{i} = kerA*P_matrices{i};
        X_projected_ker{i} = kerA'*kerA* X_matrices{i};
        X_projected_span{i} = spanA'*spanA* X_matrices{i};
        P_projected_ker{i} = kerA'*kerA* P_matrices{i};
        R_projected_span{i} = spanA'*spanA*R_matrices{i};
    end

end



%% Chybové matice
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

%Error matrices for kernel component of x_k
Error_matrices_kernel = cell(size(b,2), 1);
for k = 1:size(b,2)
    error_matrix_ker = zeros(1,maxiter);
    X = X_projected_ker{k};
    first_vector_x1 = norm(X(:,2),2);
    for i = 3:maxiter
        Norm_xi = norm(X(:,i),2);
        error_matrix_ker(1,i-2) = Norm_xi/norm(first_vector_x1,2);
        if (i > 5 && error_matrix_ker(1,i-2)==1)
            error_matrix_ker(1,i:end) =0; %pripsala jsem _ker
            break
        end
    end
    Error_matrices_kernel{k} = error_matrix_ker; %pripsala jsem _ker
end

%Error matrices for span component of x_k
Error_matrices_range = cell(size(b,2), 1);
for k = 1:size(b,2)
    error_matrix_range = zeros(1,maxiter);
    X = X_projected_span{k};
    first_vector_x1 = norm(X(:,2),2);
    for i = 1:maxiter
        Norm_xi = norm(X(:,i),2);
        error_matrix_range(1,i) = Norm_xi/norm(first_vector_x1,2);
        if (i > 5 && error_matrix_range(1,i)==1)
            error_matrix_range(1,i:end) =0;
            break
        end
    end
    Error_matrices_range{k} = error_matrix_range;
end

% checking A-orthogonality of vectors P_k against the previous vector P_{k-1}
A_orthogonality_matrices_p = cell(size(b,2), 1);
for k=1:size(b,2)
    P = P_matrices{k};
    p_a_orthogonality_matrix = zeros(1,maxiter);
    for i = 2:maxiter
        A_orthogonality_check = P(:,i)'*S*P(:,i-1);
        p_a_orthogonality_matrix(1,i-1) = abs(A_orthogonality_check);
    end
    A_orthogonality_matrices_p{k} = p_a_orthogonality_matrix;
end

%comparing Gamma components
% secist ker a span a porovnat ciste se spanem.
Span_components = cell(size(b,2), 1);
Ker_components = cell(size(b,2), 1);
Combined_components = cell(size(b,2), 1);
Denominator_values = cell(size(b,2), 1);

for k = 1: size(b,2)
    span_component_matrix = zeros(1,maxiter);
    ker_component_matrix = zeros(1,maxiter);
    combined_component_matrix = zeros(1,maxiter);
    denominator_matrix = zeros(1,maxiter);
    R = R_projected_span{k};
    P = P_matrices{k};
    b_kernel = b_ker(:,k);   
    for i=1:maxiter
        denominator = P(:,i)'*S*P(:,i);
        denominator_matrix(1,i) = denominator;
        span_component_matrix(1,i) = (R(:,i)'*R(:,i))/denominator;
        ker_component_matrix(1,i) = (b_kernel'*b_kernel)/denominator;
        combined_component_matrix(1,:) = span_component_matrix+ker_component_matrix;
    end
    Span_components{k} = span_component_matrix;
    Ker_components{k} = ker_component_matrix;
    Combined_components{k} = combined_component_matrix;
    Denominator_values{k} = denominator_matrix;
end

% checking when gamma_i/2 > gammahat
gamma_matrix_unperturbed = Gamma_matrices{1};
divergence_points = zeros(size(b,2),1);
for k = 2: size(b,2)
    gamma_matrix_perturbed = Gamma_matrices{k};
    for i = 1:maxiter
        if gamma_matrix_perturbed(1,i) > 2* gamma_matrix_unperturbed(1,i)
            divergence_points(k) = i;
            break
        else 
            if gamma_matrix_perturbed(1,i) < 0
            divergence_points(k) = i;
            break
            end
        end
    end
end


% Looking at error of line search in each step
% popremyslet, jestli to dava smysl
Error_matrices_linesearch = cell(size(b,2), 1);
for k = 1:size(b,2)
    error_matrix_linesearch = zeros(1,maxiter);
    X = X_matrices{1};
    P = P_matrices{1}; %nevim jestli mam pouzivat stejne vektory P
    gammas_linesearch = Gamma_matrices{k};
    for i = 1:maxiter
        x_linesearch = X(:,i) + gammas_linesearch(i)*P(:,i);
        A_norm_xi = sqrt((converged_x - x_linesearch)'*S*(converged_x - x_linesearch));
        error_matrix_linesearch(1,i) = A_norm_xi/(sqrt((converged_x - x0)'*S*(converged_x - x0)));
        if (i > 5 && error_matrix_linesearch(1,i)==1)
            error_matrix_linesearch(1,i:end) =0;
            break
        end
    end
    Error_matrices_linesearch{k} = error_matrix_linesearch;
end

% Looking at deviation from A-orthogonality in each step compared to
% previous direction vector from perturbed case 
A_orthogonality_deviation_matrices = cell(size(b,2), 1);
p_deviation_matrices = cell(size(b,2), 1);
for k = 1:size(b,2)
    P = P_matrices{1};
    R = R_matrices{1};
    p_a_orthogonality_deviation_matrix = zeros(1,maxiter);
    delta_deviation = Delta_matrices{k};
    p_deviation_matrices{k}(:,1) = P(:,1);
    for i = 2:maxiter
        p_deviation = R(:,i) + delta_deviation(i)*P(:,i);
        p_deviation_matrices{k}(:,i) = p_deviation; 
        % p_a_orthogonality_deviation_matrix(1,i-1) = abs((p_deviation')*S*P(:,i-1));
        p_a_orthogonality_deviation_matrix(1,i-1) = abs((p_deviation')*S*p_deviation_matrices{k}(:,i-1));
    end
    A_orthogonality_deviation_matrices{k} = p_a_orthogonality_deviation_matrix;
end


colors = {'r', 'g', 'b', 'm', 'c'};

%check dimensions
size(Error_matrices{1}(1:maxiter));
size(1:maxiter);

%% plotting

figure(3)
hold on;
for i = 1:size(b,2)
    semilogy(1:maxiter, Error_matrices{i}(1:maxiter),['o-', colors{i}]);
    hold on;
end
title('Relative error for different delta values');
xlabel('Step k');
ylabel('||x - x_i||_A / ||x - x_0||_A');
% ylim([0, inf]);
xlim([0,len]);
legend('\beta = 0', '\beta = 1e-6', '\beta = 1e-4','\beta = 1e-2', 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;

figure(4)
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

figure(5)
hold on;
for i = 1:size(b,2)
    semilogy(1:maxiter, Error_matrices_kernel{i}(1:maxiter),['o-', colors{i}]);
    hold on;
end
title('Relative error - kernel components of x_i');
xlabel('Step k');
ylabel('||x_i||/ ||x_0||');
% ylim([0, inf]);
xlim([0,len+30])
legend('\beta = 0', '\beta = 1e-6', '\beta = 1e-4','\beta = 1e-2', 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;

figure(11)
hold on;
for i = 1:size(b,2)
    semilogy(1:maxiter, Error_matrices_range{i}(1:maxiter),['o-', colors{i}]);
    hold on;
end
title('Relative error - range components of x_i');
xlabel('Step k');
ylabel('||x_i||/ ||x_0||');
% ylim([-inf, inf]);
xlim([0,len+30])
legend('\beta = 0', '\beta = 1e-6', '\beta = 1e-4','\beta = 1e-2', 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;

figure(6)
hold on;
for i = 1:size(b,2)
    semilogy(1:maxiter, A_orthogonality_matrices_p{i}(1:maxiter), ['o-', colors{i}]);
    hold on;
end
title('A orthogonality of vector p_k against previous vector p_{k-1}');
xlabel('Step k');
ylabel('A');
% ylim([0, inf]);  
xlim([0, len+30]);
legend('\beta = 0', '\beta = 1e-6', '\beta = 1e-4','\beta = 1e-2', 'Location', 'best');
set(gca, 'YScale', 'log');
yticks([-1e-13, 0, 1e-13]);
hold off;


% legend1 = ['\delta = 0', '\delta = 1e-6', '\delta = 1e-4','\delta = 1e-2', 'Location', 'best'];
% for i=1:size(b,2)
%     figure(6+i)
%     hold on;
%     semilogy(1:maxiter, Ker_components{i}(1:maxiter),'o-r');
%     hold on;
%     deltastring = legend1(i);
%     semilogy(1:maxiter, Span_components{i}(1:maxiter),'o-b');
%     title(sprintf('Comparison of Gamma components %s', deltastring));
% xlabel('Step k');
% ylabel('Value');
% ylim([0, inf]);  
% xlim([0, len]);
% set(gca, 'YScale', 'log');
% legend('Kernel Component', 'Span Component', 'Location', 'best');
% hold off;
% end


for i = 1:size(b,2)
    figure(6+i);
    hold on;
    h_combined = semilogy(1:maxiter, Combined_components{i}(1:maxiter),'.-r');
    h_span = semilogy(1:maxiter, Span_components{i}(1:maxiter),'.-b');
    title('Comparison of \gamma components');
    xlabel('Step k', 'FontName', 'AvantGarde');
    ylabel('Value', 'FontName', 'AvantGarde');
    set(gca, 'YScale', 'log');
    set([ h_combined h_span], 'LineWidth', 1.5, 'MarkerSize', 6);
    set( h_combined,'MarkerFaceColor',[.8 0 0],'MarkerEdgeColor','r');
    set( h_span, 'MarkerFaceColor',[.0 0 .8],'MarkerEdgeColor','b');
    xlim([0,len]);
    legend('Combined Components', 'Span Component', 'Location', 'best');
    hold off;
end


figure(12)
hold on;
for i = 1:size(b,2)
    semilogy(1:maxiter, Error_matrices_linesearch{i}(1:maxiter),['o-', colors{i}]);
    hold on;
end
title('Linesearch deviation in each step compared to unperturbed case.');
xlabel('Step k');
ylabel('||x - x_i||_A / ||x - x_0||_A');
% ylim([0, inf]);
xlim([0,len]);
legend('\beta = 0', '\beta = 1e-6', '\beta = 1e-4','\beta = 1e-2', 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;

figure(13)
hold on;
for i = 1:size(b,2) 
    semilogy(1:maxiter, A_orthogonality_deviation_matrices{i}(1:maxiter),['o-', colors{i}]);
    hold on;
end
title('Deviation in A-orthogonality in computation of \delta.');
xlabel('Step k');
ylabel('||x - x_i||_A / ||x - x_0||_A');
% ylim([0, inf]);
xlim([0,len]);
legend('\beta = 0', '\beta = 1e-6', '\beta = 1e-4','\beta = 1e-2', 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;