%% 1 Gradient projete
A=[1,2;-1,3];
b=[-3;4];

% 1.1 grad(J)
dJ = @(A,b,x) 2*A'*(A*x-b);
proj = @(x) max(x,zeros(size(x)));

[~,D]=eig(A'*A);
D=diag(D);
lambda_max = max(D);
lambda_min = min(D);

tau = lambda_min/lambda_max^2;
x=minimize(A,b,1.0/lambda_max-0.01,100,dJ,proj);
%% 1.7
b2=-10;
b1=10+b2;
b=[b1;b2];
x=minimize(A,b,1.0/lambda_max-0.01,100,dJ,proj);

%% 2
A=[2,1,0;1,3,1;1,0,2];
b=[3;1;3];

[V,D]=eig(A'*A);
dJ = @(A,b,x) 2*A'*(A*x-b);

D=diag(D);
lambda_max = max(D);
lambda_min = min(D);

tau = lambda_min/lambda_max^2;
x=minimize(A,b,0.089,1000,dJ,@projection);

%% 3.3
n=6;
d1 = 2*ones(n,1);
d2 = -1*ones(n-1,1);
A = diag(d1)+diag(d2,1)+diag(d2,-1);
C=[3,1,0,-1,0,0;-1,2,1,0,0,0;0,0,0,1,-1,1];
b=ones(n,1);
d=[0;1;0];

[V,D]=eig(0.5*A);
D=diag(D);
lambda_max = max(D);
lambda_min = min(D);

Csqnorm = sum(sum(C.^2));
rho = 2.1 * 2 * lambda_min / Csqnorm;
Ainv = inv(A);

lambda = zeros(3,1);
N=1000;
history = zeros(N,1);
x=zeros(n,1);
x_old = x;
for i=1:N
    x = Ainv*(b-C'*lambda);
    history(i,1) = norm(x-x_old);
    lambda = max(lambda+rho*(C*x-d),zeros(3,1));
    x_old = x;
end

plot(2:N,log(history(2:end)));
title('log |x_n-x_{n-1}|');

lambdaG = C'*lambda;
dJ = A*x - b;

%% 3.4
b=-ones(n,1);

lambda = zeros(3,1);
N=1000;
x=zeros(n,1);
x_old = x;
for i=1:N
    x = Ainv*(b-C'*lambda);
    lambda = max(lambda+rho*(C*x-d),zeros(3,1));
end

lambdaG = C'*lambda;
dJ = A*x - b;
