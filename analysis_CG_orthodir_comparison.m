% Bachelor Thesis - Conjugate Gradient Method for Solving Singular Systems
% Janika Pokorná
% May 2025, Faculty of Mathematics and Physics, Charles University

%% Run comparison of CG and Orthodir 

function[] = analysis_CG_orthodir_comparison(runName, A, b, x0, maxiter, tol, spanA, kerA)
    arguments 
        runName
        A
        b
        x0
        maxiter
        tol
        spanA
        kerA
    end

    folderPath = fullfile(pwd, runName);
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end

    b_ker = kerA'*kerA * b;


    X_matrices = cell(size(b,2), 1);
    Gamma_matrices = cell(size(b,2), 1);
    Delta_matrices = cell(size(b,2), 1);
    P_matrices = cell(size(b,2), 1);
    P_projected_ker = cell(size(b,2), 1);
    R_projected_span = cell(size(b,2), 1);
    X_projected_ker = cell(size(b,2), 1);
    R_matrices = cell(size(b,2), 1);

    for i = 1:size(b,2)
        [x, X_i,l,P_i,R_i,Gamma_i,Delta_i] = conjugate_grad(A, b(:,i),x0,maxiter,tol);
        if i == 1
            converged_x = x;
            len = l;
        end
        X_matrices{i} = X_i;
        Gamma_matrices{i} = Gamma_i;
        Delta_matrices{i} = Delta_i;
        P_matrices{i} = P_i;
        R_matrices{i} = R_i;
        
        column_norms = norm(kerA, 2);
        
        X_projected_ker{i} = kerA'*kerA* X_matrices{i};
        P_projected_ker{i} = kerA'*kerA* P_matrices{i};
        R_projected_span{i} = spanA'*spanA*R_matrices{i};
    end

    X_matrices_2 = cell(size(b,2), 1);
    Gamma_matrices_2 = cell(size(b,2), 1);
    Delta_matrices_2 = cell(size(b,2), 1);
    P_matrices_2 = cell(size(b,2), 1);
    P_projected_ker_2 = cell(size(b,2), 1);
    R_projected_span_2 = cell(size(b,2), 1);
    X_projected_ker_2 = cell(size(b,2), 1);
    R_matrices_2 = cell(size(b,2), 1);
    for i = 1:size(b,2)
    [x_2, X_i_2, l_2, P_i_2, R_i_2, Gamma_i_2, Delta_i_2] = orthodir(A, b(:,i), x0, maxiter, tol);
    if i == 1
        converged_x_2 = x_2;
        len_2 = l_2;
    end
    X_matrices_2{i} = X_i_2;
    Gamma_matrices_2{i} = Gamma_i_2;
    Delta_matrices_2{i} = Delta_i_2;
    P_matrices_2{i} = P_i_2;
    R_matrices_2{i} = R_i_2;
    column_norms = norm(kerA, 2);
    X_projected_ker_2{i} = kerA'*kerA* X_matrices_2{i};
    P_projected_ker_2{i} = kerA'*kerA* P_matrices_2{i};
    R_projected_span_2{i} = spanA'*spanA* R_matrices_2{i};
    end

    colors = {'b'};

    %% Relative error matrices
Error_matrices = cell(size(b,2), 1);
Error_matrices_2 = cell(size(b,2), 1);
for k = 1:size(b,2)
    error_matrix = zeros(1,maxiter);
    error_matrix_2 = zeros(1,maxiter);
    X = X_matrices{k};
    X_2 = X_matrices_2{k};
    for i = 1:maxiter
        A_norm_xi = sqrt((converged_x - X(:,i))'*A*(converged_x - X(:,i)));
        A_norm_xi_2 = sqrt((converged_x_2 - X_2(:,i))'*A*(converged_x_2 - X_2(:,i)));
        error_matrix(i) = A_norm_xi / sqrt((converged_x - x0)'*A*(converged_x - x0)); 
        error_matrix_2(i) = A_norm_xi_2 / sqrt((converged_x_2 - x0)'*A*(converged_x_2 - x0)); 
        if (i > 5 && error_matrix(i) == 1)
            error_matrix(i:end) = 0;
            break
        end
        if (i > 5 && error_matrix_2(i) == 1)
            error_matrix_2(i:end) = 0;
            break
        end
    end
    Error_matrices{k} = error_matrix;
    Error_matrices_2{k} = error_matrix_2;
end
    
figure
hold on;
for i = 1:size(b,2)
    % Original error curves
    h(i) = semilogy(1:maxiter, Error_matrices{i}(1:maxiter), '-', 'Color', colors{i}, 'LineWidth', 3);
    
    % _2 version curves (dashed)
    h_2(i) = semilogy(1:maxiter, Error_matrices_2{i}(1:maxiter), '--', 'Color', colors{i}, 'LineWidth', 2);
end
legend({'CG', 'Orthodir'}, 'Location', 'best');
xlim([0, len]);
set(gca, 'YScale', 'log');
set(gca, 'FontSize', 14);
hold off;

figTitle = 'Relative error for different beta values (comparison of Orthodir and CG)';
title('Relative error for different beta values (comparison of Orthodir and CG)');
safeTitle = regexprep(figTitle, '[^\w]', '_');
filename = sprintf('%s_%s.eps', runName, safeTitle);
saveas(gcf, fullfile(folderPath, filename), 'epsc');
filename = sprintf('%s_%s.png', runName, safeTitle);
saveas(gcf, fullfile(folderPath, filename));
