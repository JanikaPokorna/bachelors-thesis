data = load('matrices\regmatD2.mat', 'D');
D1 = data.D;
rank(D1)
size(D1,1)
D1 = D1'*D1;
b = D1(:,3);

[x,X] = conjugate_grad(D1,b);
error_matrix = zeros(size(X,2),1);
for i = 1:size(X,2)
    A_norm_xi = sqrt((x - X(:,i))'*D1*(x - X(:,i)));
    error_matrix(:,i) = A_norm_xi/sqrt((x - x0)'*D1*(x - x0));
end
for i = 1:size(X,2)
    semilogy(i, error_matrix(i), 'o-');
    hold on;
end
xlabel('Step k')
ylabel('||x - x_i||_A / ||x - x_0||_A')
title('Relative Error Plot')
hold off;
