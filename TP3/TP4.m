clc;
clear all;
close all;
X = double(imread('lenna.gif'));
% load gatlin2;
F = X + 10*randn(size(X));



% We test with 3 choice of lambda and 4 choices of phi
sqepilon=0;%0.01*0.01;
step=0.02;
Mu = [30,10,3,1];
M = length(Mu);

figure(2);
subplot(2,M,1)
imshow(X,[]);
title('Lena')
subplot(2,M,2)
imshow(F,[]);
title('Lena noisy')


for l = 1:M
    mu = Mu(l);
    phi = @(t) sqrt(t.^2+sqepilon);
    dphi = @(t) t./sqrt(t.^2+sqepilon);
    u1=F;
    N=300;
    history_g1 = zeros(N,1);
    history_dJ_numerical = zeros(N,1);
    history_dJ_theory = zeros(N,1);
    history_J1 = zeros(N,1);
    history_e1 = zeros(N,1);
    for i=1:N
        [gx,gy]=grad_im(u1);
        normgu = sqrt(gx.^2+gy.^2);
        dphigu=dphi(normgu);
        unit_gx = gx./(normgu+1e-7);
        unit_gy = gy./(normgu+1e-7);
        g = (u1-F)/mu - div_champ(dphigu.*unit_gx, dphigu.*unit_gy);
        
        J = sum(sum((F-u1).^2))/2/mu + sum(sum(phi(normgu)));
        u1=u1-step*g;
        history_g1(i) = norm(g);
        history_J1(i) = J;
        history_e1(i) = norm(X-u1);
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
    plot(1:N, history_e1,'DisplayName',Name);
    hold on;
    title('error u-X');
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
    subplot(2,M,l+2)
    imshow(u1,[]);
    title(Name)
end
Solu1 = u1; % mu=1

%% 2.1
close all;
mu=1;

N=200;
Tau=[0.125,0.15,0.25,0.4,0.5,1];

for t = 1:length(Tau)
    px=zeros(size(X));
    py=px;
    tau = Tau(t);
    history_dp = zeros(N,1);
    history_e = zeros(N,1);
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
    figure(1);
    subplot(121)
    plot(1:N, log(history_dp),'DisplayName',['tau=',num2str(tau)]);
    hold on;
    title('log difference p')
    legend off;
    legend show;
    
    subplot(122)
    plot(1:N, history_e,'DisplayName',['tau=',num2str(tau)]);
    hold on;
    title('log error u-X')
    legend off;
    legend show;
end

figure(2);
ProF = mu*div_champ(px,py);
u=F-ProF;
imshow([X,F,u],[])

%% 2.2
close all;
mu=1;

N=200;
Tau=[0.125,0.15,0.25];

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
        history_e(i) = norm(Solu1-u);
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
    title('log error u-solution(1)')
    legend off;
    legend show;
end

figure(2);
u=-mu*v;
imshow([X,F,u,Solu1,u-Solu1],[])
%% 3
sqsigma = 2.0;
U=-5:5;
g2 = exp(-0.5/sqsigma*U.^2)./sqrt(2*pi*sqsigma);
convA = @(u) conv2(g2,g2,u,'same');

F = convA(X) + 10*randn(size(X));


nu=0.01;
mu=1;
u=zeros(size(X));
N = 200;
M = 50;

history_e = zeros(N,1);
for i=1:N
    v = u + nu * convA(F-convA(u));
    % chambolle
    px=zeros(size(X));
    py=px;
    tau = 0.25;
    for j=1:N
        d = div_champ(px,py) - v/mu/nu;
        [gx,gy] = grad_im(d);
        denom = 1 + tau * sqrt(gx.^2+gy.^2);
        px = (px+tau*gx)./denom;
        py = (py+tau*gy)./denom;
    end
    ProF = mu*nu*div_champ(px,py);
    u=v-ProF;
    history_e(i) = norm(X-u);
end
figure(1);
plot(1:N, log(history_e));
hold on;
title('log error u-X')
legend off;
legend show;

figure(2);
imshow([X,F,u],[]);
