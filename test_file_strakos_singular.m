%% test_file strakos singular

n = 500;
ker_dim = 1;
rho = 0.8; %the smaller this is, the more eigenval are close, should be below 1
a = 5;
c = 100;
delta = 0;
[S,D,spanA,kerA] = singular_strakos(n,ker_dim,a,c,rho); % creates strakos matrix with ker of dimension 1
diagonal_values = diag(D);
% create right-hand side vector b
% b = S(:,265); %converges nicely in 30 steps
b = make_vector_b(spanA,kerA,delta);

% % Plot the diagonal values
% figure(1);
% semilogy(diagonal_values, 'or');
% title('Diagonal values of the matrix');
% xlabel('Index');
% ylabel('Diagonal Value');
% grid on;

S = S'*S;
eigen_values = eig(S);

% Plot the eigen values of symmetric PD strakos matrix 
figure(2);
semilogy(eigen_values, 'or');
title('Diagonal values of the matrix');
xlabel('Index');
ylabel('Diagonal Value');
grid on;

%check dimensions
size(b)
size(Gamma)

[x,X,len,P,R,Gamma] = conjugate_grad(S,b);

%save converged value if delta = 0;
if delta == 0
    converged_x = x;
end

%check dimensions
size(S)
size(X)
len

error_matrix = zeros(1,len);
for i = 1:len
    A_norm_xi = sqrt((converged_x - X(:,i))'*S*(converged_x - X(:,i)));
    error_matrix(1,i) = A_norm_xi/(sqrt((converged_x - x0)'*S*(converged_x - x0)));
end

%Plot of relative error
figure(3)
semilogy(1:len,error_matrix,'o-r');
xlabel('Step k')
ylabel('||x - x_i||_A / ||x - x_0||_A')
title('Relative Error Plot')
hold off;

%Plot gamma
figure(4)
hold on
plot(1:len+1,Gamma(1,1:len+1),'o-m');
xlabel('Step k')
ylabel('gamma')
title('Values of Gamma')
hold off;
