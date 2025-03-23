%% test_file strakos singular multiple plots
runName = 'Experiment_strakos';
folderPath = fullfile(pwd, runName);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

n = 500;
ker_dim = 1;
maxiter = 550;
tol = 1e-15;
x0 = zeros(n,1);
rho = 0.8; %the smaller this is, the more eigenval are close, should be below 1
a = 5;
c = 100;
betas = [0,1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2];
[S,D,spanA,kerA] = singular_strakos(n,ker_dim,a,c,rho); % creates strakos matrix with ker of dimension 1
diagonal_values = diag(D);

%% create right-hand side vector b
b = make_multi_vector_b(spanA,kerA,betas);
% b kernel component
b_ker = kerA'*kerA * b;
% b_kernel'*b_kernel;


%% Plot the eigenvalues of symmetric PD strakos matrix 
figure;
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
    Gamma_matrices{i} = Gamma_i;
    Delta_matrices{i} = Delta_i;
    P_matrices{i} = P_i;
    R_matrices{i} = R_i;
    % Oddělení komponenty v jádře
    column_norms = norm(kerA, 2);  % normy sloupců
    is_orthonormal_cols = all(abs(column_norms - 1) < eps);
    if is_orthonormal_cols == true
        %kontrola dimenzí
        size(kerA)
        size(X_matrices{i}) 
        % X_projected_ker{i} = kerA* X_matrices{i};
        % P_projected_ker{i} = kerA*P_matrices{i};
        %Ortogonálně promítnuté hodnoty:
        X_projected_ker{i} = kerA'*kerA* X_matrices{i};
        X_projected_span{i} = spanA'*spanA* X_matrices{i};
        P_projected_ker{i} = kerA'*kerA* P_matrices{i};
        R_projected_span{i} = spanA'*spanA*R_matrices{i};
    end

end



%% Matice relativní chyby aproximace
Error_matrices = cell(size(b,2), 1);
for k = 1:size(b,2)
    %error_detail_matrix = zeros(1,maxiter);
    error_matrix = zeros(1,maxiter);
    X = X_matrices{k};
    for i = 1:maxiter
        A_norm_xi = sqrt((converged_x - X(:,i))'*S*(converged_x - X(:,i)));
        error_matrix(i) = A_norm_xi/(sqrt((converged_x - x0)'*S*(converged_x - x0))); 
        if (i > 5 && error_matrix(i)==1)
            error_matrix(i:end) = 0;
            break
        end
    end
    Error_matrices{k} = error_matrix;
end

%% Vykreslení detailu divergence - je presne videt v jakém kroku metoda diverguje
Convergence_detail_matrices = cell(size(b,2), 1); %pro pozdejsi detail
for k = 1:size(b,2)
    error_detail_matrix = zeros(1,maxiter);
    error_matrix_unperturbed = Error_matrices{1};
    error_matrix_perturbed = Error_matrices{k};
    for i= 1:maxiter
        error_detail_matrix(1,i) = (error_matrix_perturbed(1,i)-error_matrix_unperturbed(1,i))/error_matrix_unperturbed(1,i);
    end
    Convergence_detail_matrices{k} = error_detail_matrix;
end

%% TOTO NENI PRILIS PRESNE, VYKRESUJEME PROJEKCI APROXIMACE (NIKOLI CHYBU)
%Error matrices for kernel component of x_k
Error_matrices_kernel = cell(size(b,2), 1);
for k = 1:size(b,2)
    error_matrix_ker = zeros(1,maxiter);
    X = X_projected_ker{k};
    first_vector_x1 = norm(X(:,2),2);
    for i = 2:maxiter
        Norm_xi = norm(X(:,i),2);
        error_matrix_ker(1,i-1) = Norm_xi/first_vector_x1;
        if (i > 5 && error_matrix_ker(1,i-1)==1)
            error_matrix_ker(1,i:end) =0;
            break
        end
    end
    Error_matrices_kernel{k} = error_matrix_ker; 
end

% TOTO NENI PRILIS PRESNE, VYKRESUJEME PROJEKCI APROXIMACE (NIKOLI CHYBU)
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
for k = 1:size(b,2)
    P = P_matrices{k};
    p_a_orthogonality_matrix = zeros(1,maxiter);
    for i = 2:maxiter
        A_orthogonality_check = P(:,i)'*S*P(:,i-1);
        p_a_orthogonality_matrix(1,i-1) = abs(A_orthogonality_check);
    end
    A_orthogonality_matrices_p{k} = p_a_orthogonality_matrix;
end

%comparing Gamma components - for comparison of Gamma/Span component of
%Gamma, and also for comparing ratio between Gamma/Gamma_hat
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
    R = R_projected_span{k};
    P = P_matrices{k};
    X = X_matrices{k};
    b_kernel = b_ker(:,k);   
    for i=1:maxiter
        denominator = P(:,i)'*S*P(:,i);
        denominator_matrix(i) = denominator;
        span_component_matrix(i) = (R(:,i)'*R(:,i))/denominator;
        ker_component_matrix(i) = (b_kernel'*b_kernel)/denominator;
        gamma_hat_matrix(i) = P(:,i)'*S*(converged_x-X(:,i))/denominator;
        combined_component_matrix(i) = span_component_matrix(i)+ker_component_matrix(i);
        %Zde se snažím spocítat Gamma/Gamma_hat:
        %OPRAVENO ZDE, PRIDAN VYPOCET i-TEHO PRVKU
        ratio_matrix(i) = combined_component_matrix(1,i)/gamma_hat_matrix(i);
        %ratio_matrix(1,i) = combined_component_matrix(1,:)/span_component_matrix(1,i);
    end
    Span_components{k} = span_component_matrix;
    Ker_components{k} = ker_component_matrix;
    Combined_components{k} = combined_component_matrix;
    Denominator_values{k} = denominator_matrix;
    Ratio_values{k} = ratio_matrix;
end

% NENI POTREBA
% checking when gamma_i/2 > gammahat - gives us a matrix of divergence
% points to check against the results of the experiment
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

% looking at ratio between gamma and gamma_hat (from unperturbed case) -
% this was the wrong way of counting the ratios
Ratio_unperturbed_values = cell(size(b,2), 1);
gamma_hat = Gamma_matrices{1};
for k = 2:size(b,2)
    gamma_per = Gamma_matrices{k};
    ratio_unperturbed_matrix = zeros(size(b,2),1);
    for i=1:maxiter
        ratio_unperturbed_matrix(1,i) = gamma_per(1,i)/gamma_hat(1,i);
    end
    Ratio_unperturbed_values{k} = ratio_unperturbed_matrix;
end


%  Looking at error of line search in each step
% nedava smysl
% Error_matrices_linesearch = cell(size(b,2), 1);
% for k = 1:size(b,2)
%     error_matrix_linesearch = zeros(1,maxiter);
%     X = X_matrices{1};
%     P = P_matrices{1}; %nevim jestli mam pouzivat stejne vektory P
%     gammas_linesearch = Gamma_matrices{k};
%     for i = 1:maxiter
%         x_linesearch = X(:,i) + gammas_linesearch(i)*P(:,i);
%         A_norm_xi = sqrt((converged_x - x_linesearch)'*S*(converged_x - x_linesearch));
%         error_matrix_linesearch(1,i) = A_norm_xi/(sqrt((converged_x - x0)'*S*(converged_x - x0)));
%         if (i > 5 && error_matrix_linesearch(1,i)==1)
%             error_matrix_linesearch(1,i:end) =0;
%             break
%         end
%     end
%     Error_matrices_linesearch{k} = error_matrix_linesearch;
% end

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
        p_deviation = R(:,i) + delta_deviation(1,i-1)*P(:,i-1);
        p_deviation_matrices{k}(:,i) = p_deviation; 
        p_a_orthogonality_deviation_matrix(1,i-1) = abs((p_deviation')*S*P(:,i-1));
        % p_a_orthogonality_deviation_matrix(1,i-1) = abs((p_deviation')*S*p_deviation_matrices{k}(:,i-1));
    end
    A_orthogonality_deviation_matrices{k} = p_a_orthogonality_deviation_matrix;
end
%%
colors = {'r', 'g', 'b', 'm', 'c', 'y', 'k', 'b'};
size(Error_matrices{1}(1:maxiter));
size(1:maxiter);

%% purely plotting values counted above
figure
hold on;
legend_entries = cell(1, size(b,2));
for i = 1:size(b,2)
    h(i) = semilogy(1:maxiter, Error_matrices{i}(1:maxiter),['o-', colors{i}]);
    hold on;
    legend_entries{i} = sprintf('\\beta = %g', betas(i));
end
title('Relative error for different delta values');
xlabel('Step k');
ylabel('||x - x_i||_A / ||x - x_0||_A');
% legend(h, legend_entries, 'Location', 'best');
xlim([0,len]);
set(gca, 'YScale', 'log');
hold off;

figTitle = get(get(gca, 'Title'), 'String');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s.epsc', runName, safeTitle);
print(gcf, '-depsc', fullfile(folderPath, filename));
%%
figure
hold on;
legend_entries = cell(1, size(b,2));
for i = 1:size(b,2)
    h(i) = plot(1:len+1, Gamma_matrices{i}(1:len+1),['o-', colors{i}]);
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
    h(i) = semilogy(1:maxiter, Error_matrices_kernel{i}(1:maxiter),['o-', colors{i}]);
    hold on;
    legend_entries{i} = sprintf('\\beta = %g', betas(i));
end
title('Relative error - kernel components of x_i');
xlabel('Step k');
ylabel('||x_i||/ ||x_0||');
% ylim([0, inf]);
xlim([0,len+30])
% legend(h, legend_entries, 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;
figTitle = get(get(gca, 'Title'), 'String');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s.epsc', runName, safeTitle);
saveas(gcf, fullfile(folderPath, filename));

figure
hold on;
legend_entries = cell(1, size(b,2));
for i = 1:size(b,2)
    h(i) = semilogy(1:maxiter, Error_matrices_range{i}(1:maxiter),['o-', colors{i}]);
    hold on;
    legend_entries{i} = sprintf('\\beta = %g', betas(i));
end
title('Relative error - range components of x_i');
xlabel('Step k');
ylabel('||x_i||/ ||x_0||');
% ylim([-inf, inf]);
xlim([0,len+30])
% legend(h, legend_entries, 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;
figTitle = get(get(gca, 'Title'), 'String');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s.epsc', runName, safeTitle);
saveas(gcf, fullfile(folderPath, filename));

figure
hold on;
legend_entries = cell(1, size(b,2));
for i = 1:size(b,2)
    h(i) = semilogy(1:maxiter, A_orthogonality_matrices_p{i}(1:maxiter), ['o-', colors{i}]);
    hold on;
    legend_entries{i} = sprintf('\\beta = %g', betas(i));
end
title('A orthogonality of vector p_k against previous vector p_{k-1}');
xlabel('Step k');
ylabel('A');
% ylim([0, inf]);  
xlim([0, len+30]);
% legend(h, legend_entries, 'Location', 'best');
set(gca, 'YScale', 'log');
yticks([-1e-13, 0, 1e-13]);
hold off;
figTitle = get(get(gca, 'Title'), 'String');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s.epsc', runName, safeTitle);
saveas(gcf, fullfile(folderPath, filename));


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
    figure
    legend_entries = cell(1, size(b,2));
    hold on;
    h_combined = semilogy(1:maxiter, Combined_components{i}(1:maxiter),'.-r');
    h_span = semilogy(1:maxiter, Span_components{i}(1:maxiter),'.-b');
    title('Comparison of gamma components');
    xlabel('Step k', 'FontName', 'AvantGarde');
    ylabel('Value', 'FontName', 'AvantGarde');
    set(gca, 'YScale', 'log');
    set([ h_combined h_span], 'LineWidth', 1.5, 'MarkerSize', 6);
    set( h_combined,'MarkerFaceColor',[.8 0 0],'MarkerEdgeColor','r');
    set( h_span, 'MarkerFaceColor',[.0 0 .8],'MarkerEdgeColor','b');
    xlim([0,len]);
    legend('Combined Components', 'Span Component', 'Location', 'best');
    hold off;
    figTitle = get(get(gca, 'Title'), 'String');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s_%d.epsc', runName, safeTitle, i);
saveas(gcf, fullfile(folderPath, filename));
end

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

figure;
hold on;
legend_entries = cell(1, size(b,2));
for i = 1:size(b,2) 
    semilogy(1:maxiter, A_orthogonality_deviation_matrices{i}(1:maxiter),['o-', colors{i}]);
    hold on;
    legend_entries{i} = sprintf('\\beta = %g', betas(i));
end
title('Deviation in A-orthogonality in computation of \delta.');
xlabel('Step k');
ylabel('||x - x_i||_A / ||x - x_0||_A');
% ylim([0, inf]);
xlim([0,len]);
% legend(h, legend_entries, 'Location', 'best');
set(gca, 'YScale', 'log');
hold off;
figTitle = get(get(gca, 'Title'), 'String');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s.epsc', runName, safeTitle);
saveas(gcf, fullfile(folderPath, filename));
%%
figure 
for i = 1:size(b,2)
    hold on;
    h_ratio = plot(1:maxiter,Convergence_detail_matrices{i}(1:maxiter),['o-', colors{i}]);
    legend('Ratio of a-norm of x in perturbed vs. unperturbed ', 'Location', 'best');
    xlim([0,len]);
    ylim([0,1.5]);
    hold off;
end
title('The Ratio of convergence.');
figTitle = get(get(gca, 'Title'), 'String');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s.epsc', runName, safeTitle);
saveas(gcf, fullfile(folderPath, filename));