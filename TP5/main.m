%% 1 Please run each section one by one without cleaning data
% Leman FENG flm8620@gmail.com
clc;
clear all;
close all;
X = double(imread('cameraman.tif'));
v = X + 10*randn(size(X));

Lambda = 10.^(-3:0.1:3);
M = length(Lambda);
figure(1);
subplot(1,3,1)
imshow(X,[]);
title('im')
subplot(1,3,2)
imshow(v,[]);
title('im noisy')
camrestb10 = cell(M);
errors = zeros(M,1);
Niter=30;
tau=0.25;

h = waitbar(0,'Initializing waitbar...');
for l = 1:M
    waitbar(l/M,h,sprintf('%d/%d', l, M));
    u=Chambolle(v,Lambda(l),Niter,tau);
    camrestb10{l}=u;
    errors(l)=sum(sum((X-u).^2));
    if l==1
        best = u;
        best_error = errors(l);
    elseif best_error > errors(l)
        best_error = errors(l);
        best = u;
    end
end
delete(h);
figure(2)
plot(log10(Lambda),errors);
xlabel('log10 lambda');
ylabel('sq error')

figure(1)
subplot(1,3,3)
imshow(best,[]);
title('best')

% the error function decreases then increases
[~,l_best] = min(errors);
best_lambda = Lambda(l_best);
fprintf('best lambda = %f\n',best_lambda)

%% 2.1
close all;
% calculate O(lambda, sigma_0, N_0) table
sigma_0 = 10;
N_0 = 256*256;
O = zeros(M,1);
b=sigma_0*randn(256, 256);
v=b;
h = waitbar(0,'Initializing waitbar...');
for l = 1:M
    waitbar(l/M,h,sprintf('%d/%d', l, M));
    u=Chambolle(v,Lambda(l),Niter,tau);
    O(l)=2 * sum(sum(u.*b));
end
delete(h);
figure(3)
plot(log10(Lambda),O);
hold on;
xlabel('log10 lambda');
ylabel('O(lambda)')
title('O(lambda, sigma0, N0) table')
% O(l,s,kN) = k*O(l,s,N) simply because the inner product is taken on a
% larger grid
% O(kl,ks,N) = k*k*O(l,s,N) because:
% minimizing F_{lambda}(u,v) = |u-v|^2 + lambda*TV(u) is equivalent to
% minimizing k^2*F_{lambda}(u,v) = |ku-kv|^2 + k*lambda*TV(ku)
% k^2*F_{lambda}(u) can be seen as F_{k*lambda}(ku,kv)
% where kv has k times noise than v, and since u_optim is solution for 
% min_{u} F_{lambda}(u,v), so k*u_optim is solution for 
% min_{u} F_{k*lambda}(u,kv)
% so <u_optim, b> is k^2 times larger. Proved.

%% 2.2
close all;
imnames={'lenna.gif','cameraman.tif','rice.png'};
for sigma = [5,10,15]
    for k=1:length(imnames)
        Otest = zeros(M,1);
        u0 = double(imread(imnames{k}));
        b = sigma*randn(size(X));
        v = u0+b;
        N = numel(v);
        
        h = waitbar(0,'Initializing waitbar...');
        for l = 1:M
            waitbar(l/M,h,sprintf('%d/%d', l, M));
            u=Chambolle(v,Lambda(l),Niter,tau);
            Otest(l)=2 * sum(sum(u.*b));
        end
        delete(h);
        % I plot Otest table for each test image and rescale it to compare later to
        % standard O table
        figure(4)
        plot(log10(Lambda*sigma_0/sigma),Otest*sigma_0^2*N_0/sigma^2/N,'DisplayName',[imnames{k},' sigma=',num2str(sigma)]);
        xlabel('log10 lambda');
        ylabel('O(lambda)')
        hold on;
        legend off;
        legend show;
    end
end
% Compare to standard table O(lambda, sigma_0, N_0)
figure(4)
plot(log10(Lambda),O,'-.', 'DisplayName','standard');
xlabel('log10 lambda');
ylabel('O(lambda)')
title('Compare to standard O(lambda, sigma0, N0) table')
hold on;
legend off;
legend show;
% It seems function O only differs a constant for different images and
% sigma, but when lambda is large, standard O table evaluated from constant
% images seems different from O calculated from real images.


%% 3
close all;
sigma = 10;
u0 = double(imread('cameraman.tif'));
b=sigma*randn(size(X));
v = u0 + b;
N = numel(v);

errors = zeros(M,1);
sq_u_lambda_v = zeros(M,1);
O_estimate = zeros(M,1);
O_real = zeros(M,1);

h = waitbar(0,'Initializing waitbar...');
for l = 1:M
    waitbar(l/M,h,sprintf('%d/%d', l, M));
    u=Chambolle(v,Lambda(l)*sigma/sigma_0,Niter,tau);
    O_estimate(l) = sigma^2*N/sigma_0^2/N_0 * O(l);
    errors(l)=sum(sum((X-u).^2));
    sq_u_lambda_v(l) = sum(sum((v-u).^2));
end
delete(h);

vb2 = 2 * sum(sum(v.*b));
bb = sum(sum(b.*b));
figure(5)
plot(log10(Lambda*sigma/sigma_0), errors, 'DisplayName', 'true error')
hold on
plot(log10(Lambda*sigma/sigma_0), O_estimate+sq_u_lambda_v, 'DisplayName', 'Estimated O + |u-v|^2')
legend off;
legend show;

[~,l_estimated_best] = min(O_estimate+sq_u_lambda_v);
estimated_best_lambda = Lambda(l_estimated_best)*sigma/sigma_0;
[~,l_true_best] = min(errors);
true_best_lambda = Lambda(l_true_best)*sigma/sigma_0;
fprintf('true best lambda = %f, estimated lambda = %f\n',true_best_lambda, estimated_best_lambda)

figure(6)
u=Chambolle(v,estimated_best_lambda,Niter,tau);
imshow([u0,v,u],[])

%% 3.2
close all;
sigma = 10;
u0 = double(imread('cameraman.tif'));
b=sigma*randn(size(X));
v = u0 + b;
N = numel(v);

% five points bisection method to find a minimum of a function which
% decrease then increase

evaluate_error = @(l) sum(sum((v-Chambolle(v,Lambda(l)*sigma/sigma_0,Niter,tau) ).^2)) + sigma^2*N/sigma_0^2/N_0 * O(l);

[argmin, minvalue, tried_x, tried_f] = five_point_bisection_minimum(evaluate_error, 1, M);

figure(7)
scatter(log10(Lambda(tried_x)*sigma/sigma_0), tried_f);
title('evaluated values in five points bisection method')
best_lambda_found =  Lambda(argmin)*sigma/sigma_0;
fprintf('best lambda found = %f\n',best_lambda_found)

figure(8)
u=Chambolle(v,estimated_best_lambda,Niter,tau);
imshow([u0,v,u],[])

