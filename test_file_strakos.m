%% test_file strakos n = 500
% this test went well, relative error converges 'linearly' in semilogy and
% shows nice curve in regular plot. converges in 89 steps.
n = 500;
rho = 0.8; %the smaller this is, the more eigenval are close, should be below 1
a = 1;
c = 100;
[S,D] = strakos(n,a,c,rho);
diagonal_values = diag(D); %eigenvalues of strakos matrix
b = rand(n,1);

% Plot the diagonal values
figure(1);
semilogy(diagonal_values, 'or');
title('Diagonal values of the matrix');
xlabel('Index');
ylabel('Diagonal Value');
grid on;

S = S'*S;
eigen_values = eig(S);

% Plot the eigen values of symmetric PD strakos matrix 
figure(2);
semilogy(eigen_values, 'or');
title('Diagonal values of the matrix');
xlabel('Index');
ylabel('Diagonal Value');
grid on;

[x,X,len] = conjugate_grad(S,b);
x;
error_matrix = zeros(1,len);
for i = 1:len 
    A_norm_xi = sqrt((x - X(:,i))'*S*(x - X(:,i)));
    error_matrix(1,i) = A_norm_xi/(sqrt((x - x0)'*S*(x - x0)));
end
%Plot of relative error
figure(3)
plot(1:len,error_matrix,'o-r');
xlabel('Step k')
ylabel('||x - x_i||_A / ||x - x_0||_A')
title('Relative Error Plot')
hold off;

