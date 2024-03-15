% data = load('matrices\regmatD2.mat', 'D');
% D1 = data.D;
% rank(D1)
% size(D1,1)
% D1 = D1'*D1;
% b = D1(:,3);
% 
% [x,X,len] = conjugate_grad(D1,b);
% x;
% error_matrix = zeros(1,len);
% for i = 1:len
%     A_norm_xi = sqrt((x - X(:,i))'*D1*(x - X(:,i)));
%     error_matrix(1,i) = A_norm_xi/sqrt((x - x0)'*D1*(x - x0));
% end
% semilogy(1:len,error_matrix,'o-r');
% xlabel('Step k')
% ylabel('||x - x_i||_A / ||x - x_0||_A')
% title('Relative Error Plot')
% hold off;

%%
% data = load('matrices\regmatD1.mat', 'D');
% D1 = data.D;
% rank(D1)
% size(D1,1)
% D1 = D1'*D1;
% b = D1(:,3);
% 
% [x,X,len] = conjugate_grad(D1,b);
% x;
% error_matrix = zeros(1,len);
% for i = 1:len
%     A_norm_xi = sqrt((x - X(:,i))'*D1*(x - X(:,i)));
%     error_matrix(1,i) = A_norm_xi/sqrt((x - x0)'*D1*(x - x0));
% end
% semilogy(1:len,error_matrix,'o-r');
% xlabel('Step k')
% ylabel('||x - x_i||_A / ||x - x_0||_A')
% title('Relative Error Plot')
% hold off;
%%
% data = load('matrices\regmatD3.mat', 'D');
% D1 = data.D;
% rank(D1)
% size(D1,1)
% D1 = D1'*D1;
% b = ones(size(D1,1),1);
% x0 = zeros(size(D1,1),1);
% 
% [x,X,len] = conjugate_grad(D1,b);
% x
% error_matrix = zeros(1,len);
% for i = 1:len
%     A_norm_xi = sqrt((x - X(:,i))'*D1*(x - X(:,i)));
%     error_matrix(1,i) = A_norm_xi/(sqrt((x - x0)'*D1*(x - x0)));
% end
% semilogy(1:len,error_matrix,'o-r');
% xlabel('Step k')
% ylabel('||x - x_i||_A / ||x - x_0||_A')
% title('Relative Error Plot')
% hold off;
%%
% data = load('matrices\bcsstm06.mat', 'Problem'); %JAK NACIST
% D1 = data.Problem.A;
% b = ones(size(D1,1),1);
% x0 = zeros(size(D1,1),1);
% 
% [x,X,len] = conjugate_grad(D1,b);
% x;
% error_matrix = zeros(1,len);
% for i = 1:len
%     A_norm_xi = sqrt((x - X(:,i))'*D1*(x - X(:,i)));
%     error_matrix(1,i) = A_norm_xi/(sqrt((x - x0)'*D1*(x - x0)));
% end
% semilogy(1:len,error_matrix,'o-r');
% xlabel('Step k')
% ylabel('||x - x_i||_A / ||x - x_0||_A')
% title('Relative Error Plot')
% hold off;

%%
% data = load('matrices\singmatS1.mat', 'S');
% D1 = data.S;
% rank(D1)
% size(D1,1)
% D1 = D1'*D1;
% b = D1(:,3);
% x0 = zeros(size(D1,1),1);
% 
% [x,X,len] = conjugate_grad(D1,b);
% x;
% error_matrix = zeros(1,len);
% size(error_matrix)
% for i = 1:len %% incompatible array sizes
%     A_norm_xi = sqrt((x - X(:,i))'*D1*(x - X(:,i)));
%     error_matrix(1,i) = A_norm_xi/(sqrt((x - x0)'*D1*(x - x0)));
% end
% semilogy(1:len,error_matrix,'o-r');
% xlabel('Step k')
% ylabel('||x - x_i||_A / ||x - x_0||_A')
% title('Relative Error Plot')
% hold off;
%%
% data = load('matrices\singmatS2.mat', 'S');
% D1 = data.S;
% rank(D1)
% size(D1,1)
% D1 = D1'*D1;
% b = D1(:,3);
% x0 = zeros(size(D1,1),1);
% 
% [x,X,len] = conjugate_grad(D1,b);
% x;
% error_matrix = zeros(1,len);
% for i = 1:len 
%     A_norm_xi = sqrt((x - X(:,i))'*D1*(x - X(:,i)));
%     error_matrix(1,i) = A_norm_xi/(sqrt((x - x0)'*D1*(x - x0)));
% end
% semilogy(1:len,error_matrix,'o-r');
% xlabel('Step k')
% ylabel('||x - x_i||_A / ||x - x_0||_A')
% title('Relative Error Plot')
% hold off;
%%
data = load('singmatS6.mat', 'S');
spanA = load('spanS6.mat','spanA') ;
kerA = load('kerS6.mat','kerA');
spanA = spanA.spanA;
kerA = kerA.kerA;
D1 = data.S;
% rank(D1)
% size(D1,1)
D1 = D1'*D1;

eigen_values = eig(D1)

% Plot the diagonal values
figure(1);
stem(eigen_values, 'filled');
title('Diagonal values of the matrix');
xlabel('Index');
ylabel('Diagonal Value');
grid on;


b = make_vector_b(spanA,kerA,1e-100);
x0 = zeros(size(D1,1),1);

[x,X,len] = conjugate_grad(D1,b); %% warning, imaginary parts of complex x or y ignored
x;
error_matrix = zeros(1,len);
size(error_matrix);
for i = 1:len 
    A_norm_xi = sqrt((x - X(:,i))'*D1*(x - X(:,i)));
    error_matrix(1,i) = A_norm_xi/(sqrt((x - x0)'*D1*(x - x0)));
end
figure(2);
semilogy(1:len,error_matrix,'o-r');
xlabel('Step k')
ylabel('||x - x_i||_A / ||x - x_0||_A')
title('Relative Error Plot')