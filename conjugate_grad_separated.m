%% CG separated into spanA/kerA

function [x,X_s,X_k,i,P_s,P_k,R_s,R_k,Gamma] = conjugate_grad_separated(A, b, spanA, kerA, x0, maxiter, tol)
    arguments
        A
        b
        spanA
        kerA
        x0 = zeros(size(A, 1),1)
        maxiter = 100;
        tol = 1e-7
    end
    n = size(A, 1); 
    R_s = zeros(n,maxiter+1);
    R_k = zeros(n,maxiter+1);
    P_s = zeros(n, maxiter+1);
    P_k = zeros(n, maxiter+1);
    X_s = zeros(n, maxiter+1);
    X_k = zeros(n, maxiter+1);
    Gamma = zeros(1,maxiter+1);
    % projecting onto spanA and kerA
    S = spanA * pinv(spanA);
    K = kerA * pinv(kerA);
    R_s(:,1) = S * (b - A * x0);
    R_k(:,1) = K * (b - A * x0);
    P_s(:,1) = R_s(:,1);
    P_k(:,1) = R_k(:,1);
    X_s(:,1) = S * x0;
    X_k(:,1) = K * x0;
    r_prod_old = (R_s(:,1)+R_k(:,1))' * (R_s(:,1)+R_k(:,1)); %%%% tady jsem skoncila zatim
    i = 1;
    while (sqrt(r_prod_old) > tol && i < maxiter)
        r_prod_old = (R_s(:,i)+R_k(:,i))' * (R_s(:,i)+R_k(:,i));
        p_prod =  A * (P_s(:,i)+P_k(:,i));
        Gamma(1,i) = r_prod_old/((P_s(:,i)+P_k(:,i))' * p_prod);
        R_s(:,i+1) = R_s(:,i) - Gamma(1,i) * p_prod;
        R_k(:,i+1) = R_k(:,i) - Gamma(1,i) * p_prod;
        X_s(:,i+1) = X_s(:,i) + Gamma(1,i) * P_s(:,i);
        X_k(:,i+1) = X_k(:,i) + Gamma(1,i) * P_k(:,i);
        r_prod_new = (R_s(:,i+1)+R_k(:,i+1))' * (R_s(:,i+1)+R_k(:,i+1));
        delta = r_prod_new / r_prod_old;
        P_s(:,i+1) = R_s(:,i+1) + delta * P_s(:,i);
        P_k(:,i+1) = R_k(:,i+1) + delta * P_k(:,i);
        r_prod_old = r_prod_new;
        i = i+1;
    end
    x = X_s(:,i) +  X_k(:,i);
    if i < maxiter
        text = sprintf('Method converged in %d iterations.\n', i);
        disp(text);
    else
        text = sprintf('Method was stopped after %d iterations.\n', i);
        disp(text);
    end
    
end