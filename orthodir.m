%% Orthodir

function [x,X,i,P,R,Gamma,Delta] = orthodir(A, b, x0, maxiter, tol)
    arguments
        A
        b
        x0 = zeros(size(A, 1),1)
        maxiter = 100;
        tol = 1e-7
    end
    n = size(A, 1); 
    S = zeros(n,maxiter+1);
    mu = zeros(1,maxiter+1);
    sigma = zeros(1,maxiter+1);
    R = zeros(n,maxiter+1);
    P = zeros(n, maxiter+2);
    X = zeros(n, maxiter+1);
    Gamma = zeros(1,maxiter+1);
    Delta = zeros(1,maxiter+1);
    R(:,2) = b;
    P(:,2) = A*b;
    x = x0;
    r_prod_old = R(:,2)' * R(:,2);
    i = 2;
    while (sqrt(r_prod_old) > tol && i < maxiter)
        r_p_prod = R(:,i)' * P(:,i);
        S(:,i) =  A * P(:,i);
        mu(1,i) = sqrt(P(:,i)' * S(:,i));
        Gamma(1,i) = r_p_prod/(P(:,i)' * S(:,i));
        R(:,i+1) = R(:,i) - Gamma(1,i) * S(:,i);
        X(:,i) = x;
        x = x + Gamma(1,i) * P(:,i);
        r_prod_new = S(:,i)' * S(:,i);
        Delta(1,i) = r_prod_new / (P(:,i)' * S(:,i));
        if i==2
            sigma(1,i) = 0;
        else 
            sigma(1,i) = mu(1,i)/(mu(1,i)-1);
        end
        P(:,i+1) = S(:,i)/mu(1,i) - (Delta(1,i)/mu(1,i)) *P(:,i) - sigma(1,i) * P(:,i-1);
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


