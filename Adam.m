close all
clear all
%% read the data
A = readmatrix("Adjacency_matrix.csv");
%% set up initial parameters
N=size(A,1);
tol=10^(-11);
% number of maximum iteration
kmax=10^5;
%% initial positions of the vertices
x=normrnd(0,N,[N,1]);
y=normrnd(0,N,[N,1]);
w=[x;y];

%% optimization using Adam
lam = 0.01;
espi = 1e-8;
m = 0;
v = 0;
k = 1;
beta1 = 0.9;
beta2 = 0.999;
while k<=kmax
    alpha = 5/(1+0.01*k);
    g = -forces(x,y,A);
    gnorm_list(k) = norm(g);
    if mod(k,200)==0
        fprintf('at iteration %d, gnorm: |g| = %.4e\n',k,norm(g));
    end
    if norm(g)<tol
        break
    end

    m = beta1*m + (1-beta1)*g;
    v = beta2*v + (1-beta2)*(g.*g);
    mhat = m/(1-beta1^k);
    vhat = v/(1-beta2^k);
    w = w - alpha*mhat./(sqrt(vhat)+espi);
    x=w(1:N);
    y=w(N+1:2*N);
    k = k+1;
end

%% plot the graph
x=w(1:N);
y=w(N+1:2*N);
plot_graph(x,y,A);
figure(2);
plot(log(gnorm_list));
xlabel('iteration');
ylabel('norm of the force in log scale');




