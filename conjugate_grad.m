% Bachelor Thesis - Conjugate Gradient Method for Solving Singular Systems
% Janika Pokorná
% May 2025, Faculty of Mathematics and Physics, Charles University

%% Conjugate Gradient method

function [x,X,i,P,R,Gamma,Delta] = conjugate_grad(A, b, x0, maxiter, tol)
    arguments
        A
        b
        x0 = zeros(size(A, 1),1)
        maxiter = 100;
        tol = 1e-7
    end
    n = size(A, 1); 
    R = zeros(n,maxiter+1);
    P = zeros(n, maxiter+1);
    X = zeros(n, maxiter+1);
    Gamma = zeros(1,maxiter+1);
    Delta = zeros(1,maxiter+1);
    R(:,1) = b - A * x0;
    P(:,1) = R(:,1);
    x = x0;
    r_prod_old = R(:,1)' * R(:,1);
    i = 1;
    while (sqrt(r_prod_old) > tol && i < maxiter)
        r_prod_old = R(:,i)' * R(:,i);
        p_prod =  A * P(:,i);
        Gamma(1,i) = r_prod_old/(P(:,i)' * p_prod);
        R(:,i+1) = R(:,i) - Gamma(1,i) * p_prod;
        X(:,i) = x;
        x = x + Gamma(1,i) * P(:,i);
        r_prod_new = R(:,i+1)' * R(:,i+1);
        Delta(1,i) = r_prod_new / r_prod_old;
        P(:,i+1) = R(:,i+1) + Delta(1,i) * P(:,i);
        r_prod_old = r_prod_new;
        i = i+1;
    end
    X(:,i) = x;
    if i < maxiter
        text = sprintf('Method converged in %d iterations.\n', i);
        disp(text);
    else
        text = sprintf('Method was stopped after %d iterations.\n', i);
        disp(text);
    end
    
end



