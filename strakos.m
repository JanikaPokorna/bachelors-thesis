function [A,D] = strakos(n,a,b,rho)

%    
    lambda      = zeros(n,1);
    lambda(1)   = a;
    lambda(n)   = b;

%
    for i=2:n-1
        lambda(i)=a+((i-1)/(n-1))*(b-a)*rho^(n-i);
    end
%    
    D = diag(lambda);
    R = rand(n);
    [Q,~] = qr(R);
    A = Q'*D*Q;