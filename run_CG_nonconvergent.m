function[] = run_CG_nonconvergent(runName, A, b, x0, maxiter, tol, betas, spanA, kerA)
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
            A_norm_xi = sqrt((X(:,i))'*A*( X(:,i)));
            error_matrix(i) = A_norm_xi/(sqrt((x0)'*A*(x0))); 
            if (i > 5 && error_matrix(i)==1)
                error_matrix(i:end) = 0;
                break
            end
        end
        Error_matrices{k} = error_matrix;
    end

%kontrola abs hodnoty p_k^TAp_k
    Checking_range_p = cell(size(b,2), 1);
    for k = 1:size(b,2)
        P = P_matrices{k};
        p_A_p_value_matrix = zeros(maxiter);
        for i = 1:maxiter
            p_A_p_value_matrix(i) = abs(P(:,i)'*A*P(:,i));
        end
        Checking_range_p{k} = p_A_p_value_matrix;
    end

 %% Vykresleni
 colors = {'r', 'g', 'b', 'm', 'c', 'k', 'b', 'b'};
 figure
    hold on;
    legend_entries = cell(1, size(b,2));
    for i = 1:size(b,2)
        h(i) = semilogy(1:maxiter, Error_matrices{i}(1:maxiter),['o-', colors{i}]);
        hold on;
        xlim([0, len]);
        legend_entries{i} = sprintf('\\beta = %g', betas(i));
    end
    legend( legend_entries, 'Location', 'best');
    xlim([0,30]);
    set(gca, 'YScale', 'log');
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
        h(i) = semilogy(1:maxiter, Checking_range_p{i}(1:maxiter), ['-', colors{i}]);
        hold on;
        xlim([0, 60]);
        legend_entries{i} = sprintf('\\beta = %g', betas(i));
    end
    
    legend(legend_entries, 'Location', 'best');
    set(gca, 'YScale', 'log');
    hold off;
    figTitle = 'Deviation in A-orthogonality in computation of \delta.';
    safeTitle = regexprep(figTitle, '[^\w]', '_');
    filename = sprintf('%s_%s.eps', runName, safeTitle);
    saveas(gcf, fullfile(folderPath, filename),'epsc');
    filename = sprintf('%s_%s.png', runName, safeTitle);
    saveas(gcf, fullfile(folderPath, filename));
end
