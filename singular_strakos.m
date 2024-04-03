function [A,D,spanA, kerA] = singular_strakos(n, kernel_dim, a, b, rho)
    % 
    lambda = zeros(n, 1);
    lambda(1) = a;
    lambda(n + 1 - kernel_dim:end) = 0; 

    % 
    for i = 2:(n - kernel_dim)
        lambda(i) = a + ((i-1)/(n-1)) * (b-a) * rho^(n-i+1);
    end

    % 
    D = diag(lambda);
    R = rand(n);
    [Q,~] = qr(R);
    spanA = Q(1:end-kernel_dim, :);
    kerA = Q(end-kernel_dim+1:end, :);
    A = Q'*D*Q;
end