%% PCG

%  Let us have a system of linear equations Ax = b for vector x, A is a
%  positive-definite Hermitian n x n matrix, and a preconditioning matrix
%  C.

function [x] = precond_conjugate_grad(A, b, C, x0, maxiter, tol)
    arguments
        A
        b
        C = eye(size(A, 1))
        x0 = zeros(size(A, 1),1)
        maxiter = 1000
        tol = 1e-4
    end
    n = size(A, 1); 
    R_ = zeros(n,maxiter+1);
    P_ = zeros(n, maxiter+1);
    R_(:,1) = b - A * x0;
    y = C\R_(:,1);
    Z(:,1) = C'\y;
    P_(:,1) = Z(:,1);
    x_ = x0;
    r_prod_old = Z(:,1)' * R_(:,1);
    iteration = 0;
    i = 1;
    while (sqrt(r_prod_old) > tol && iteration < maxiter)
        r_prod_old = Z(:,i)' * R_(:,i);
        p_prod =  A * P_(:,i);
        gamma = r_prod_old/(P_(:,i)' * p_prod);
        x_ = x_ + gamma * P_(:,i);
        R_(:,i+1) = R_(:,i) - gamma * p_prod;
        y = C\R_(:,i+1);
        Z(:,i+1) = C'\y;
        r_prod_new = Z(:,i+1)' * R_(:,i+1);
        delta = r_prod_new / r_prod_old;
        P_(:,i+1) = Z(:,i+1) + delta * P_(:,i);
        r_prod_old = r_prod_new;
        iteration = iteration + 1;
        i = i+1;
    end
    R = R_;
    P = P_;
    x = x_;
    if iteration < maxiter
        text = sprintf('Method converged in %d iterations.\n', iteration);
        disp(text);
    else
        text = sprintf('Method was stopped after %d iterations.\n', iteration);
        disp(text);
    end
end
