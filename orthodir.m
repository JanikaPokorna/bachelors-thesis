%% Orthodir

function [x,X,i,P,R,Gamma,Delta] = orthodir(A, b, x0, maxiter, tol)
    arguments
        A
        b
        x0
        maxiter = 100;
        tol = 1e-4
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
    i = 1;
    R(:,1) = b;
    P(:,1) = A*b;
    while (norm(R(:,i)) > tol && i < maxiter+2)
        S(:,i) =  A * P(:,i);
        mu(i) = sqrt(P(:,i)' * S(:,i));
        Gamma(i) =  ( R(:,i)' * P(:,i) )/(mu(i)^2);
        X(:,i+1) = X(:,i) + Gamma(i) * P(:,i);
        R(:,i+1) = R(:,i) - Gamma(i) * S(:,i);
         if i==1
            sigma(i)=0;
            third_term = 0;
        else 
            sigma(i) = mu(i)/(mu(i-1));
            third_term = sigma(i)*P(:,i-1);
        end
        Delta(i) =  ( S(:,i)' * S(:,i) )/ (mu(i)^2);
        P(:,i+1) = S(:,i)/mu(i) - (Delta(i)/mu(i)) *P(:,i) - third_term;
        i = i+1;
    end
    
    if i < maxiter
        text = sprintf('Method converged in %d iterations.\n', i);
        disp(text);
    else
        text = sprintf('Method was stopped after %d iterations.\n', i);
        disp(text);
    end
    x = X(:,i);
    i = i-2;
end
