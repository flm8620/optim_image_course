%% 3.1
clc;
clear all;
close all;

%load gatlin2;
X = double(imread('lenna.gif'));
% X = X/max(max(X));
[gx,gy]=grad_im(X);
figure(1);
imshow(X,[]);
title('im');
colorbar;

figure(2);
imshow([gx,gy],[]);
title('grad');

d = div_champ(gx,gy);
figure(3);
imshow([X,d],[]);
title('im & laplacien');

%F = imnoise(X,'gaussian',0,0.001);
F = X + 10*randn(size(X));

gradF = @(f,u,l) 2*(u-f)-2*l*laplacien_im(u);
u=zeros(size(F));
lambda = 2.5;
N = 250;
history_g = zeros(N,1);
history_J = zeros(N,1);
for i=1:N
    g = gradF(F,u,lambda);
    u=u-0.01*g;
    [gx,gy]=grad_im(u);
    J = sum(sum((F-u).^2)) + lambda * sum(sum(gx.^2+gy.^2));
    history_g(i) = norm(g);
    history_J(i) = J;
end
figure(5);
imshow([X,F,u],[]);
title('left to right: original, noised, recovered')
figure(6);
plot(1:N, log(history_g));
title('log gradient norm')
figure(7);
plot(1:N, log(history_J));
title('log energy')
%% 3.2
clc;
clear all;
close all;
rng(123);
X = double(imread('lenna.gif'));
F = X + 10*randn(size(X));
S = fsym(F);
[M,N]=size(S);
lambda = 2.5;

FS = fft2(S);
[p,q]=meshgrid(0:N-1, 0:M-1);
Fu = FS./(1+4*lambda*(sin(pi*p/N).^2 + sin(pi*q/M).^2));
u = ifft2(Fu);
u = real(unsym(u));

figure(10);
imshow([X,F,u],[]);
title('original, noisy, and recovered images');

