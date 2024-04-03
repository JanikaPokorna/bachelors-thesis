n=100;
rho = 2;
a = 0.1;
b = 100;
[S,D] = strakos(n,a,b,rho);
diagonal_values = diag(D);

% Plot the diagonal values
figure;
stem(diagonal_values, 'filled');
title('Diagonal values of matrix S');
xlabel('Index');
ylabel('Diagonal Value');
grid on;