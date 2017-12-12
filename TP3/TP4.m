%% 1
clc;
clear all;
close all;
X = double(imread('lenna.gif'));
% load gatlin2;
F = X + 10*randn(size(X));

% We test with 3 choice of mu
sqepilon=0.0000;%0.01*0.01;
step=0.02;
Mu = [1,3,6,10];
M = length(Mu);

figure(2);
subplot(2,M,2)
imshow(X,[]);
title('Lena')
subplot(2,M,3)
imshow(F,[]);
title('Lena noisy')

h = waitbar(0,'Initializing waitbar...');
for l = 1:M
    mu = Mu(l);
    phi = @(t) sqrt(t.^2+sqepilon);
    dphi = @(t) t./sqrt(t.^2+sqepilon);
    u=F;
    N=400;
    history_g1 = zeros(N,1);
    history_dJ_numerical = zeros(N,1);
    history_dJ_theory = zeros(N,1);
    history_J1 = zeros(N,1);
    history_e1 = zeros(N,1);
    waitbar(l/M,h,sprintf('%d/%d', l, M));
    for i=1:N
        [gx,gy]=grad_im(u);
        normgu = sqrt(gx.^2+gy.^2);
        dphigu=dphi(normgu);
        unit_gx = gx./(normgu+1e-7);
        unit_gy = gy./(normgu+1e-7);
        g = (u-F)/mu - div_champ(dphigu.*unit_gx, dphigu.*unit_gy);
        
        J = sum(sum((F-u).^2))/2/mu + sum(sum(phi(normgu)));
        u=u-step*g;
        history_g1(i) = norm(g);
        history_J1(i) = J;
        history_e1(i) = norm(X-u);
        if i>1
            history_dJ_theory(i-1)=-step*sum(sum(g.*g));
            history_dJ_numerical(i-1)=history_J1(i)-history_J1(i-1);
        end
    end
    
    Name=['mu = ',num2str(mu)];
    figure(1);
    subplot(141)
    plot(1:N, log(history_g1),'DisplayName',Name);
    hold on;
    title('log gradient norm')
    legend off;
    legend show;
    subplot(142)
    plot(1:N, log(history_J1),'DisplayName',Name);
    hold on;
    title('log energy')
    legend off;
    legend show;
    subplot(143)
    plot(1:N, log(history_e1),'DisplayName',Name);
    hold on;
    title('log error u-X');
    legend off;
    legend show;
    subplot(144)
    plot(1:N, history_dJ_numerical,'DisplayName',[Name,' numer'],'LineStyle','-');
    hold on;
    plot(1:N, history_dJ_theory,'DisplayName',[Name,' theo'],'LineStyle','--');
    title('dJ theory & numerical');
    legend off;
    legend show;
    
    figure(2);
    subplot(2,M,l+M)
    imshow(u,[]);
    title(Name)
end
delete(h);
Solu1 = u; % mu=10


%% 2.1
close all;
mu=10;

N=400;
Tau=[20,5,0.5,0.25,0.125];
M = length(Tau);
figure(2);
subplot(2,M,2)
imshow(X,[]);
title('Lena')
subplot(2,M,3)
imshow(F,[]);
title('Lena noisy')

h = waitbar(0,'Initializing waitbar...');
for t = 1:length(Tau)
    px=zeros(size(X));
    py=px;
    tau = Tau(t);
    history_dp = zeros(N,1);
    history_e = zeros(N,1);
    waitbar(t/M,h,sprintf('%d/%d', t, M));
    for i=1:N
        d = div_champ(px,py) - F/mu;
        [gx,gy] = grad_im(d);
        denom = 1 + tau * sqrt(gx.^2+gy.^2);
        
        px1 = (px+tau*gx)./denom;
        py1 = (py+tau*gy)./denom;
        history_dp(i)=sum(sum((px1-px).^2+(py1-py).^2));
        ProF = mu*div_champ(px,py);
        u=F-ProF;
        history_e(i) = norm(X-u);
        px=px1;
        py=py1;
    end
    Name = ['tau=',num2str(tau)];
    figure(1);
    subplot(121)
    plot(1:N, log(history_dp),'DisplayName',Name);
    hold on;
    title('log difference p')
    legend off;
    legend show;
    
    subplot(122)
    plot(1:N, log(history_e),'DisplayName',Name);
    hold on;
    title('log error u-X')
    legend off;
    legend show;
    
    ProF = mu*div_champ(px,py);
    u=F-ProF;
    figure(2);
    subplot(2,M,t+M)
    imshow(u,[]);
    title(Name)
end

delete(h);

figure(3);
ProF = mu*div_champ(px,py);
Solu21=F-ProF;
subplot(2,3,1)
imshow(X,[]);
colorbar;
title('Lena')
subplot(2,3,2)
imshow(F,[]);
title('Lena noisy')
subplot(2,3,4)
imshow(Solu21,[]);
title('Chambolle')
subplot(2,3,5)
imshow(Solu1,[]);
title('Gradient descent')
subplot(2,3,6)
imshow(Solu21-Solu1,[])
colorbar;
title('difference')

%% 2.2
close all;
mu=10;

N=400;
Tau=[0.5,0.3,0.25,0.15,0.125];

for t = 1:length(Tau)
    px=zeros(size(X));
    py=px;
    tau = Tau(t);
    history_dp = zeros(N,1);
    history_e = zeros(N,1);
    for i=1:N
        v = -F/mu + div_champ(px,py);
        [gx,gy] = grad_im(v);
        
        px1 = (px+tau*gx);
        py1 = (py+tau*gy);
        
        denom = max(1, sqrt(px1.^2+py1.^2));
        px1=px1./denom;
        py1=py1./denom;
        
        history_dp(i)=sum(sum((px1-px).^2+(py1-py).^2));
        u=-mu*v;
        history_e(i) = norm(Solu21-u);
        px=px1;
        py=py1;
    end
    figure(1);
    subplot(121)
    plot(1:N, log(history_dp),'DisplayName',['tau=',num2str(tau)]);
    hold on;
    title('log difference p')
    legend off;
    legend show;
    
    subplot(122)
    plot(1:N, log(history_e),'DisplayName',['tau=',num2str(tau)]);
    hold on;
    title('log error solution(2.2)-solution(2.1)')
    legend off;
    legend show;
end

Solu22=-mu*v;
figure(2);

subplot(2,3,1)
imshow(X,[]);colorbar;
title('Lena')
subplot(2,3,2)
imshow(F,[]);
title('Lena noisy')
subplot(2,3,4)
imshow(Solu21,[]);
title('Chambolle')
subplot(2,3,5)
imshow(Solu22,[]);
title('Projected Gradient');
subplot(2,3,6)
imshow(Solu22-Solu21,[])
colorbar;
title('difference is negligible')
%% 3
close all;
sqsigma = 2.0;
U=-5:5;
g2 = exp(-0.5/sqsigma*U.^2)./sqrt(2*pi*sqsigma);
convA = @(u) conv2(g2,g2,u,'same');

F = convA(X) + 10*randn(size(X));
%F = X + 10*randn(size(X));


nu=0.1;
mu=10;
u=F;
N = 200;
M = 5;

history_e = zeros(N*M,1);
h = waitbar(0,'Initializing waitbar...');
for i=1:N
    waitbar(i/N,h,sprintf('%d/%d', i, N));
    v = u + nu * convA(F-convA(u));
    % chambolle
    px=zeros(size(X));
    py=px;
    tau = 0.25;
    for j=1:M
        d = div_champ(px,py) - v/mu/nu;
        [gx,gy] = grad_im(d);
        denom = 1 + tau * sqrt(gx.^2+gy.^2);
        px = (px+tau*gx)./denom;
        py = (py+tau*gy)./denom;
        ProF = mu*nu*div_champ(px,py);
        u=v-ProF;
        history_e((i-1)*M+j) = norm(X-u);
    end
end
delete(h);
Decon = u;

% gradient descent
% we use the code in Question 4 in last TP+
u=F;
h = waitbar(0,'Initializing waitbar...');
history_e1 = zeros(N*M,1);
step=0.005;
for i=1:N*M
    waitbar(i/N/M,h,sprintf('%d/%d', i, N*M));
    % deconvolution
    [gx,gy]=grad_im(u);
    normgu = sqrt(gx.^2+gy.^2);
    % gradient
    g1 = convA(convA(u)-F)/mu - div_champ(gx, gy);
    u=u-step*g1;
    history_e1(i) = norm(X-u);
end


delete(h);
gradient_descent = u;

figure(1)
plot(1:N*M, log(history_e1),'DisplayName','Gradient descent');
hold on;
plot(1:N*M, log(history_e),'DisplayName','Deconvolution');
title('log error u-X')
legend off;
legend show;


figure(2);
subplot(1,4,1)
imshow(X,[0,255]);
title('Lena')
subplot(1,4,2)
imshow(F,[0,255]);
title('Lena noisy')
subplot(1,4,3)
imshow(Decon,[0,255]);
title('Deconvolution')
subplot(1,4,4)
imshow(gradient_descent,[0,255]);
title('Gradient descent')
