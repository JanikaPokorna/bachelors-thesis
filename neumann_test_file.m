%% Neumann test file
runName = 'Experiment_neumann';

folderPath = fullfile(pwd, runName);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

%parameters
n = 25;
sizeA = n^2;
A = gallery('neumann',sizeA);
maxiter = 1000;
tol = 1e-15;
x0 = zeros(sizeA,1);
betas = [0,1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2];

%right hand side vector b
[V,D] = eigs(A,sizeA);
span_A = V(:,1:sizeA-1);
ker_A = V(:,sizeA);
b = make_multi_vector_b(span_A',ker_A',betas); % musim do funkce dosadit radky
b_ker = ker_A*ker_A' * b;

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
Delta_matrices = cell(size(b,2), 1);
P_matrices = cell(size(b,2), 1);
P_projected_ker = cell(size(b,2), 1);
X_projected_ker = cell(size(b,2), 1);
R_matrices = cell(size(b,2), 1);
R_projected_span = cell(size(b,2), 1);
X_projected_span = cell(size(b,2), 1);

for i = 1:size(b,2)
    [x, X_i,l,P_i,R_i,Gamma_i,Delta_i] = conjugate_grad(A, b(:,i),x0,maxiter,tol);
    if i == 1
        converged_x = x;
        len = l;
    end
    X_matrices{i} = X_i;
    Gamma_matrices{i,1} = Gamma_i;
    P_matrices{i,1} = P_i;
    R_matrices{i,1} = R_i;
    % Oddělení komponenty v jádře
    column_norms = norm(ker_A', 2);  % normy sloupců
    is_orthonormal_cols = all(abs(column_norms - 1) < eps);
    if is_orthonormal_cols == true
        %kontrola dimenzí
        size(ker_A')
        size(X_matrices{i}) 
    end
    X_projected_ker{i} = (ker_A*ker_A')*X_matrices{i};
    P_projected_ker{i} = (ker_A*ker_A')*P_matrices{i};
    R_projected_span{i} = (span_A*span_A')*R_matrices{i};
end

%% error matrices

Error_matrices = cell(size(b,2), 1);
for k = 1:size(b,2)
    error_matrix = zeros(1,maxiter);
    X = X_matrices{k};
    for i = 1:maxiter
        M_norm_xi = sqrt((converged_x - X(:,i))'*A*(converged_x - X(:,i)));
        error_matrix(1,i) = M_norm_xi/(sqrt((converged_x - x0)'*A*(converged_x - x0)));
        if (i > 5 && error_matrix(1,i)==1)
            error_matrix(1,i:end) = 0;
            break
        end
    end
    Error_matrices{k} = error_matrix;
end

%% Vykreslení detailu divergence
Convergence_detail_matrices = cell(size(b,2), 1); %pro pozdejsi detail
for k=1:size(b,2)
    error_detail_matrix = zeros(1,maxiter);
    error_matrix_unperturbed = Error_matrices{1};
    error_matrix_perturbed = Error_matrices{k};
    for i= 1:maxiter
        error_detail_matrix(1,i) = (error_matrix_perturbed(1,i)-error_matrix_unperturbed(1,i))/error_matrix_unperturbed(1,i);
    end
    Convergence_detail_matrices{k} = error_detail_matrix;
end

%Error matrices for kernel component of x_k
% nedava smysl vykreslovat relativne, protoze delime skoro nulou
Error_matrices_kernel = cell(size(b,2), 1);
for k = 1:size(b,2)
    error_matrix_ker = zeros(1,maxiter);
    X = X_projected_ker{k};
    first_vector_x1 = norm(X(:,2),2);
    for i = 2:maxiter
        Norm_xi = norm(X(:,i),2);
        error_matrix_ker(1,i-1) = Norm_xi;
        % error_matrix_ker(1,i-1) = Norm_xi/first_vector_x1;
        if (i > 5 && error_matrix_ker(1,i-1)==1)
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
Combined_components = cell(size(b,2), 1);
Denominator_values = cell(size(b,2), 1);
Ratio_values = cell(size(b,2), 1);

for k = 1: size(b,2)
    span_component_matrix = zeros(1,maxiter);
    ker_component_matrix = zeros(1,maxiter);
    combined_component_matrix = zeros(1,maxiter);
    denominator_matrix = zeros(1,maxiter);
    gamma_hat_matrix = zeros(1,maxiter);
    ratio_matrix = zeros(1,maxiter);
    R = R_matrices{k};
    P = P_matrices{k};
    b_kernel = b_ker(:,k);
    for i=1:maxiter
        denominator = P(:,i)'*A*P(:,i);
        span_component_matrix(1,i) = (R(:,i)'*R(:,i))/denominator;
        ker_component_matrix(1,i) = (b_kernel'*b_kernel)/denominator;
        gamma_hat_matrix(i) = P(:,i)'*A*(converged_x-X(:,i))/denominator;
        combined_component_matrix(i) = span_component_matrix(i)+ker_component_matrix(i);
        ratio_matrix(i) = combined_component_matrix(1,i)/gamma_hat_matrix(i);
    end
    Span_components{k} = span_component_matrix;
    Ker_components{k} = ker_component_matrix;
    Combined_components{k} = combined_component_matrix;
    Denominator_values{k} = denominator_matrix;
    Ratio_values{k} = ratio_matrix;
end

%%
divergence_points_approximate = zeros(1, size(b,2));
for k = 1:size(b,2)
    difference = 0;
    error_matrix = Error_matrices{k};
    error_matrix_unperturbed = Error_matrices{1};
    for i=1:maxiter
        if (i>5 && error_matrix(1,i) > 1.59*error_matrix_unperturbed(1,i))
            divergence_points_approximate(1,k) = i;
            break
        end
    end
end

%Searching for max value of gamma/gamma_hat ratio before divergence
max_ratio_values = zeros(1, size(b,2));
 for k = 1:size(b,2)
    divergence_point = divergence_points_approximate(1,k);
    lower_bound = max(1, divergence_point - 500);
    upper_bound = min(divergence_point+500,maxiter); 
    max_ratio_values(1,k) = max(Ratio_values{k}(lower_bound:upper_bound));
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
        end
    end
end

colors = {'r', 'g', 'b', 'm', 'c', 'y', 'k', 'g', 'b'};

%check dimensions
size(Error_matrices{1}(1:maxiter));
size(1:maxiter);

%% plotting 

figure
hold on;
legend_entries = cell(1, size(b,2));
for i = 1:size(b,2)
    h(i) = semilogy(1:maxiter, Error_matrices{i}(1:maxiter),['o-', colors{i}]);
    hold on;
    legend_entries{i} = sprintf('\\beta = %g', betas(i));
end
title('Relative error for different beta values');
xlabel('Step k');
ylabel('||x - x_i||_A / ||x - x_0||_A');
ylim([0, inf]);
xlim([0,len]);
legend(legend_entries, 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;
figTitle = get(get(gca, 'Title'), 'String');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s.epsc', runName, safeTitle);
saveas(gcf, fullfile(folderPath, filename));
%%
figure
hold on;
legend_entries = cell(1, size(b,2));
for i = 1:size(b,2)
    h(i) =semilogy(1:len+1, Gamma_matrices{i}(1:len+1),['o-', colors{i}]);
    hold on;
    legend_entries{i} = sprintf('\\beta = %g', betas(i));
end
xlabel('Step k')
ylabel('Gamma values')
title('Values of Gamma')
% legend(h, legend_entries, 'Location', 'best');
hold off;
figTitle = get(get(gca, 'Title'), 'String');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s.epsc', runName, safeTitle);
saveas(gcf, fullfile(folderPath, filename));

figure
hold on;
legend_entries = cell(1, size(b,2));
for i = 1:size(b,2)
    h(i)= semilogy(1:maxiter, Error_matrices_kernel{i}(1:maxiter),['o-', colors{i}]);
    hold on;
    legend_entries{i} = sprintf('\\beta = %g', betas(i));
end
title('Relative error - kernel components of x_i');
xlabel('Step k');
ylabel('||x_i||/ ||x_0||');
ylim([0, inf]);
xlim([0,len+30])
% legend(h, legend_entries, 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;
figTitle = get(get(gca, 'Title'), 'String');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s.epsc', runName, safeTitle);
saveas(gcf, fullfile(folderPath, filename));

% for i = 1:size(b,2)
%     figure;
%     hold on;
%     h_combined = semilogy(1:maxiter, Combined_components{i}(1:maxiter),'.-r');
%     h_span = semilogy(1:maxiter, Span_components{i}(1:maxiter),'.-b');
%     title('Comparison of \gamma components');
%     xlabel('Step k', 'FontName', 'AvantGarde');
%     ylabel('Value', 'FontName', 'AvantGarde');
%     set(gca, 'YScale', 'log');
%     set([ h_combined h_span], 'LineWidth', 1.5, 'MarkerSize', 6);
%     set( h_combined,'MarkerFaceColor',[.8 0 0],'MarkerEdgeColor','r');
%     set( h_span, 'MarkerFaceColor',[.0 0 .8],'MarkerEdgeColor','b');
%     xlim([divergence_points_approximate(1,i)-100,divergence_points_approximate(1,i)]);
%     legend('Combined Components', 'Span Component', 'Location', 'best');
%     hold off;
% end

for i = 1:size(b,2)
    figure
    hold on;
    h_relative = plot(1:maxiter,Ratio_values{i}(1:maxiter),'.-b');
    plot(1:maxiter,2,'.-r');
    legend('Ratio of gamma and gamma_hat', 'Location', 'best');
    xlim([0,len]);
    ylim([0,3]);
    hold off;
    figTitle = get(get(gca, 'Title'), 'String');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s_%d.epsc', runName, safeTitle, i);
saveas(gcf, fullfile(folderPath, filename));
end

figure 
for i = 1:size(b,2)
    hold on;
    h_ratio = plot(1:maxiter,Convergence_detail_matrices{i}(1:maxiter),'.-b');
    legend('Ratio of a-norm of x in perturbed vs. unperturbed ', 'Location', 'best');
    xlim([0,len]);
    ylim([0,1.5]);
    hold off;
end
figTitle = get(get(gca, 'Title'), 'String');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s.epsc', runName, safeTitle);
saveas(gcf, fullfile(folderPath, filename));