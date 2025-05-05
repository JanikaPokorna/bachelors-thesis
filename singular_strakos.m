% Bachelor Thesis - Conjugate Gradient Method for Solving Singular Systems
% Janika Pokorná
% May 2025, Faculty of Mathematics and Physics, Charles University

function [A,D,spanA, kerA] = singular_strakos(n, kernel_dim, a, b, rho)
    % 
    lambda = zeros(n, 1);
    lambda(1) = a;
    if kernel_dim ~= 0
    lambda(n + 1 - kernel_dim:end) = 0; 
    end
    

    % 
    for i = 2:(n - kernel_dim)
        lambda(i) = a + ((i-1)/(n-1)) * (b-a) * rho^(n-i+1);
    end

    % 
    D = diag(lambda);
    
    figure;
    semilogy(diag(D), 'or');
    grid on;

    R = rand(n);
    [Q,~] = qr(R);
    spanA = Q(1:end-kernel_dim, :);
    if kernel_dim ~= 0
        kerA = Q(end-kernel_dim+1:end, :); %matice s dim = kernle_dim x n, plus už je ortogonální (tzn sloupce jsou ortonormální)
    else 
        kerA = zeros(1,n);
    end
    
    A = Q'*D*Q; %uz je sym
end