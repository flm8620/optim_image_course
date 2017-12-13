clc;
clear all;
close all;
X = double(imread('cameraman.tif'));
F = X + 10*randn(size(X));

Lambda = 10.^(-3:0.1:3);
M = length(Lambda);
figure(1);
subplot(1,3,1)
imshow(X,[]);
title('im')
subplot(1,3,2)
imshow(F,[]);
title('im noisy')
camrestb10 = cell(M);
errors = zeros(M,1);
Niter=30;
tau=0.25;

h = waitbar(0,'Initializing waitbar...');
for l = 1:M
    waitbar(l/M,h,sprintf('%d/%d', l, M));
    u=Chambolle(F,Lambda(l),Niter,tau);
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

%% 2.1
% calculate O(lambda, sigma_0, N_0) table
sigma_0 = 10;
N_0 = 256*256;
O = zeros(M,1);
imnames={'lenna.gif','cameraman.tif','rice.png'};
for k=1:length(imnames)
    for a=1:3
        X = double(imread(imnames{k}));
        u0= zeros(size(X));
        b=sigma_0*randn(size(X));
        v=u0+b;
        
        h = waitbar(0,'Initializing waitbar...');
        for l = 1:M
            waitbar(l/M,h,sprintf('%d/%d', l, M));
            u=Chambolle(v,Lambda(l),Niter,tau);
            O(l)=2 * sum(sum(u.*b));
        end
        delete(h);
        
        figure(3)
        plot(log10(Lambda),log(O),'DisplayName',[imnames{k},' trial ', num2str(a)]);
        xlabel('log10 lambda');
        ylabel('log O(lambda)')
        hold on;
        legend off;
        legend show;
    end
end
figure(3)
title('O(lambda,sigma0,N0) for different images and different noise')

%% 2.2
Otest = zeros(M,1);
sigma = 5;
imnames={'circuit.tif'};
for k=1:length(imnames)
    X = double(imread(imnames{k}));
    u0= zeros(size(X));
    b=sigma*randn(size(X));
    v=u0+b;
    N = numel(X);
    
    h = waitbar(0,'Initializing waitbar...');
    for l = 1:M
        waitbar(l/M,h,sprintf('%d/%d', l, M));
        u=Chambolle(v,Lambda(l),Niter,tau);
        Otest(l)=2 * sum(sum(u.*b));
    end
    delete(h);
    
    figure(3)
    plot(log10(Lambda),Otest,'DisplayName',imnames{k});
    xlabel('log10 lambda');
    ylabel('O(lambda)')
    hold on;
    legend off;
    legend show;
end

% Compare to O(lambda, sigma_0, N_0)
figure(3)
plot(log10(Lambda*sigma/sigma_0),O*sigma^2*N/sigma_0^2/N_0,'DisplayName','standard');
xlabel('log10 lambda');
ylabel('O(lambda)')
hold on;
legend off;
legend show;

%% 3


