function[] = run_analysis_Orthodir(runName, A, b, x0, maxiter, tol, betas, spanA, kerA)
    arguments 
        runName
        A
        b
        x0
        maxiter
        tol
        betas
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
    %CG + zaznamenání do proměnných
    for i = 1:size(b,2)
        [x, X_i,l,P_i,R_i,Gamma_i,Delta_i] = orthodir(A, b(:,i),x0,maxiter,tol);
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
            
        end
        X_projected_ker{i} = kerA'*kerA* X_matrices{i};
        P_projected_ker{i} = kerA'*kerA* P_matrices{i};
        R_projected_span{i} = spanA'*spanA*R_matrices{i};
    end
    
    
    
    %% Matice relativní chyby aproximace
    Error_matrices = cell(size(b,2), 1);
    for k = 1:size(b,2)
        error_matrix = zeros(1,maxiter);
        X = X_matrices{k};
        for i = 1:maxiter
            A_norm_xi = sqrt((converged_x - X(:,i))'*A*(converged_x - X(:,i)));
            error_matrix(i) = A_norm_xi/(sqrt((converged_x - x0)'*A*(converged_x - x0))); 
            if (i > 5 && error_matrix(i)==1)
                error_matrix(i:end) = 0;
                break
            end
        end
        Error_matrices{k} = error_matrix;
    end
    
    %% Vykreslení detailu divergence
    Divergence_detail_matrices = cell(size(b,2), 1);
    for k = 1:size(b,2)
        error_detail_matrix = zeros(1,maxiter);
        error_matrix_unperturbed = Error_matrices{1};
        error_matrix_perturbed = Error_matrices{k};
        for i= 1:maxiter
            error_detail_matrix(1,i) = (error_matrix_perturbed(1,i)-error_matrix_unperturbed(1,i))/error_matrix_unperturbed(1,i);
        end
        Divergence_detail_matrices{k} = error_detail_matrix;
    end
    
    %Error matrices for kernel component of x_k
    Error_matrices_kernel = cell(size(b,2), 1);
    for k = 1:size(b,2)
        error_matrix_ker = zeros(1,maxiter);
        X = X_projected_ker{k};
        for i = 2:maxiter
            error_matrix_ker(1,i-1) = norm(X(:,i),2);
            if (i > 5 && error_matrix_ker(1,i-1)==1)
                error_matrix_ker(1,i:end-1) = 0;
                break
            end
        end
        Error_matrices_kernel{k} = error_matrix_ker; 
    end
    



    % % TOTO NENI PRILIS PRESNE, VYKRESUJEME PROJEKCI APROXIMACE (NIKOLI CHYBU)
    % %Error matrices for span component of x_k
    % Error_matrices_range = cell(size(b,2), 1);
    % for k = 1:size(b,2)
    %     error_matrix_range = zeros(1,maxiter);
    %     X = X_projected_span{k};
    %     for i = 1:maxiter
    %         Norm_xi = norm(X(:,i),2);
    %         error_matrix_range(1,i) = Norm_xi;
    %         if (i > 5 && error_matrix_range(1,i)==1)
    %             error_matrix_range(1,i:end) =0;
    %             break
    %         end
    %     end
    %     Error_matrices_range{k} = error_matrix_range;
    % end
    
    % checking A-orthogonality of vectors P_k against the previous vector P_{k-1}
    A_orthogonality_matrices_p = cell(size(b,2), 1);
    for k = 1:size(b,2)
        P = P_matrices{k};
        p_a_orthogonality_matrix = zeros(1,maxiter);
        for i = 2:maxiter
            A_orthogonality_check = P(:,i)'*A*P(:,i-1);
            p_a_orthogonality_matrix(1,i-1) = abs(A_orthogonality_check);
        end
        A_orthogonality_matrices_p{k} = p_a_orthogonality_matrix;
    end
    
   
    colors = {'r', 'g', 'b', 'm', 'c', 'k', 'b', 'b'};
    size(Error_matrices{1}(1:maxiter));
    size(1:maxiter);
    
    %% purely plotting values counted above
    figure
    hold on;
    legend_entries = cell(1, size(b,2));
    for i = 1:size(b,2)
        h(i) = semilogy(1:maxiter, Error_matrices{i}(1:maxiter),'-', 'Color', colors{i}, 'LineWidth', 3);
        hold on;
        legend_entries{i} = sprintf('\\beta = %g', betas(i));
    end
    legend( legend_entries, 'Location', 'best');
    xlim([0,len]);
    set(gca, 'YScale', 'log');
    set(gca, 'FontSize', 14);
    hold off;
    figTitle = 'Relative error for different beta values';
    safeTitle = regexprep(figTitle, '[^\w]', '_');
    filename = sprintf('%s_%s.eps', runName, safeTitle);
    saveas(gcf, fullfile(folderPath, filename),'epsc');
    filename = sprintf('%s_%s.png', runName, safeTitle);
    saveas(gcf, fullfile(folderPath, filename));
   
    
    figure
    hold on;
    legend_entries = cell(1, size(b,2));
    for i = 1:size(b,2)
        h(i) = semilogy(1:maxiter, Error_matrices_kernel{i}(1:maxiter),'-', 'Color', colors{i}, 'LineWidth', 3);
        hold on;
        legend_entries{i} = sprintf('\\beta = %g', betas(i));
    end
    % ylim([0, inf]);
    xlim([0,len]);
    legend(legend_entries, 'Location', 'best');
    set(gca, 'YScale', 'log');
    set(gca, 'FontSize', 14);
    hold off;
    figTitle = 'Relative error - kernel components of x_i';
    safeTitle = regexprep(figTitle, '[^\w]', '_');
    filename = sprintf('%s_%s.eps', runName, safeTitle);
    saveas(gcf, fullfile(folderPath, filename),'epsc');
    filename = sprintf('%s_%s.png', runName, safeTitle);
    saveas(gcf, fullfile(folderPath, filename));
   
    figure
    hold on;
    legend_entries = cell(1, size(b,2));
    for i = 1:size(b,2)
        h(i) = semilogy(1:maxiter, A_orthogonality_matrices_p{i}(1:maxiter), '-', 'Color', colors{i}, 'LineWidth', 3);
        hold on;
        legend_entries{i} = sprintf('\\beta = %g', betas(i));
    end
    xlim([0, len]);
    legend(legend_entries, 'Location', 'best');
    set(gca, 'YScale', 'log');
    set(gca, 'FontSize', 14);
    hold off;
    figTitle = 'Deviation in A-orthogonality in computation of \delta.';
    safeTitle = regexprep(figTitle, '[^\w]', '_');
    filename = sprintf('%s_%s.eps', runName, safeTitle);
    saveas(gcf, fullfile(folderPath, filename),'epsc');
    filename = sprintf('%s_%s.png', runName, safeTitle);
    saveas(gcf, fullfile(folderPath, filename));
  
    
    figure 
    for i = 1:size(b,2)
        hold on;
        h_ratio = plot(1:maxiter,Divergence_detail_matrices{i}(1:maxiter),'-', 'Color', colors{i}, 'LineWidth', 3);
        xlim([0,len]);
        set(gca, 'FontSize', 14);
        ylim([0,1.5]);
        hold off;
    end
    figTitle = 'The Detail of divergence.';
    safeTitle = regexprep(figTitle, '[^\w]', '_');
    filename = sprintf('%s_%s.eps', runName, safeTitle);
    saveas(gcf, fullfile(folderPath, filename),'epsc');
    filename = sprintf('%s_%s.png', runName, safeTitle);
    saveas(gcf, fullfile(folderPath, filename));
end
