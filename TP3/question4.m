%%
% X = double(imread('lenna.gif'));
load gatlin2;
F = X + 3*randn(size(X));


%phi = @(t) t;
%dphi = @(t) ones(size(t));
phi = @(t) t.^2./(1+t.^2);
dphi = @(t) (2*t.*(1+t.^2)-2*t.^3)./(1+t.^2).^2;
% phi = @(t) t.^2;
% dphi = @(t) 2*t;

% phi = @(t) log(1+t.^2);
% dphi = @(t) 2*t./(1+t.^2);


phi = @(t) 2*sqrt(1+t.^2)-2;
dphi = @(t) 2*t./sqrt(1+t.^2);

gradF = @(f,u,l) 2*(u-f)-2*l*div_champ(dphi()*u);
u=F;
lambda = 1;
N=1000;
history_g = zeros(N,1);
history_J = zeros(N,1);
for i=1:N
    [gx,gy]=grad_im(u);
    normgu = sqrt(gx.^2+gy.^2);
    dphigu=dphi(normgu);
    unit_gx = gx./(normgu+1e-7);
    unit_gy = gy./(normgu+1e-7);
    g = 2*(u-F) - lambda*div_champ(dphigu.*unit_gx, dphigu.*unit_gy);

    J = sum(sum((F-u).^2)) + lambda * sum(sum(phi(normgu)));
    u=u-0.01*g;
    history_g(i) = norm(g);
    history_J(i) = J;
end

figure(1);
plot(1:N, log(history_g));
title('gradient norm')
figure(2);
plot(1:N, log(history_J));
title('energy')
figure(3);
imshow([X,F,u],[]);
