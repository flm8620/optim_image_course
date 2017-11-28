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
figure(4);
imshow(F,[]);
colorbar;

gradF = @(f,u,l) 2*(u-f)-2*l*laplacien_im(u);
u=zeros(size(F));
lambda = 5;
N = 1000;
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
imshow([F,u],[]);

figure(6);
plot(1:N, history_g);
title('gradient norm')
figure(7);
plot(1:N, history_J);
title('energy')
%%
rng(123);
load gatlin2;
X = double(imread('lenna.gif'));
%X = double(sum(imread('len_top.jpg'),3))/3;
F = X + 10*randn(size(X));
S = fsym(F);
[M,N]=size(S);
% imshow(S,[]);
lambda = 10;

FS = fft2(S);

figure(8);
imshow(fftshift(angle(FS)),[-pi,pi]);
colormap(gca,hsv);
colorbar;
title('angle of Fourier of sym image')

%[p,q]=meshgrid([0:N/2 -N/2+1:-1], [0:M/2 -M/2+1:-1]);
[p,q]=meshgrid(0:N-1, 0:M-1);
Fu = FS./(1+4*lambda*(sin(pi*p/N).^2 + sin(pi*q/M).^2));
u = ifft2(Fu);


figure(9);
imshow(angle(u),[-pi,pi]);
colormap(gca,hsv);
colorbar;
title('angle of u')

u = real(unsym(u));
figure(10);
imshow([X,F,u],[]);
title('original, noisy, and recovered images');

