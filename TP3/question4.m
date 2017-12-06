%% Question 4.1
clc;
clear all;
close all;
% X = double(imread('lenna.gif'));
load gatlin2;
F = X + 3*randn(size(X));

Phi = {@(t) t;
    @(t) t.^2./(1+t.^2);
    @(t) log(1+t.^2);
    @(t) 2*sqrt(1+t.^2)-2;
    };
dPhi = {@(t) ones(size(t));
    @(t) (2*t.*(1+t.^2)-2*t.^3)./(1+t.^2).^2;
    @(t) 2*t./(1+t.^2);
    @(t) 2*t./sqrt(1+t.^2);
    };
Names={'t','t^2/(1+t^2)','log(1+t^2)','2sqrt(1+t^2)-2'};

% We test with 3 choice of lambda and 4 choices of phi
rows={[],[],[]};
linS = {'-','--',':'};
Lambda = [10,3,1];
for l = 1:3
    lambda = Lambda(l);
    u={[],[],[],[]};
    for k=1:4
        phi = Phi{k};
        dphi = dPhi{k};
        u1=F;
        N=100;
        history_g1 = zeros(N,1);
        history_J1 = zeros(N,1);
        history_e1 = zeros(N,1);
        for i=1:N
            [gx,gy]=grad_im(u1);
            normgu = sqrt(gx.^2+gy.^2);
            dphigu=dphi(normgu);
            unit_gx = gx./(normgu+1e-7);
            unit_gy = gy./(normgu+1e-7);
            g = 2*(u1-F) - lambda*div_champ(dphigu.*unit_gx, dphigu.*unit_gy);
            
            J = sum(sum((F-u1).^2)) + lambda * sum(sum(phi(normgu)));
            u1=u1-0.03*g;
            history_g1(i) = norm(g);
            history_J1(i) = J;
            history_e1(i) = norm(X-u1);
        end
        u{k}=u1;
        figure(1);
        subplot(131)
        plot(1:N, log(history_g1),'DisplayName',Names{k}, 'LineStyle',linS{l});
        hold on;
        title('log gradient norm')
        legend show;
        subplot(132)
        plot(1:N, log(history_J1),'DisplayName',Names{k}, 'LineStyle',linS{l});
        hold on;
        title('log energy')
        legend show;
        subplot(133)
        plot(1:N, history_e1,'DisplayName',Names{k}, 'LineStyle',linS{l});
        hold on;
        title('error');
        legend show;
        
    end
    rows{l} = [u{1},u{2},u{3},u{4}];
end
figure(3);
imshow([rows{1};rows{2};rows{3}],[]);
title('top to buttom: lambda=10,3,1   left to right: 4 choices of phi')
% comment:
% it seems choices with 
% {log(1+t^2),     lambda=10}
% {t,              lambda=3}
% {2sqrt(1+t^2)-2, lambda = 1}
% have better performance. we should adapt choice of lambdaof with phi

% the concativity of phi is important.
% convex phi is more like a t^2 function
% linear or concave phi is more like total variance minimization

%% Question 4.2
clc;
clear all;
close all;
% It's easy to deduce that the adjoint operator for a convolution with
% kernal g(x,y) is the convolution with kernel h(x,y):=g(x,y).
% Since gaussian kernel is symmetric, so gaussian convolution is
% self-adjoint operator.
% So the derivative of ||f-A(u)||^2 w.r.t u is 2A(A(u) - f)

% now we use lenna as a test image X, then we convolve X with gaussian, at
% last we add noise.

X = double(imread('lenna.gif'));
sqsigma = 2.0;
U=-5:5;
g2 = exp(-0.5/sqsigma*U.^2)./sqrt(2*pi*sqsigma);
convA = @(u) conv2(g2,g2,u,'same');

F = convA(X);
F = F + 5*randn(size(X));

% We use phi = t, total variation minimization
phi = @(t) t;
dphi = @(t) ones(size(t));

% We compare Denoising and Deconvolution
% u1 for deconvolution
% u2 for denoising
u1=F;
u2=F;
lambda = 3; % change lambda and run again this section to have multiple plot
N=200;
history_g1 = zeros(N,1);
history_J1 = zeros(N,1);
history_e1 = zeros(N,1); % error too original image
history_g2 = zeros(N,1);
history_J2 = zeros(N,1);
history_e2 = zeros(N,1);

step=0.05;
for i=1:N
    % deconvolution
    [gx,gy]=grad_im(u1);
    normgu = sqrt(gx.^2+gy.^2);
    dphigu=dphi(normgu);
    unit_gx = gx./(normgu+1e-7);
    unit_gy = gy./(normgu+1e-7);
    % gradient
    g1 = 2*convA(convA(u1)-F) - lambda*div_champ(dphigu.*unit_gx, dphigu.*unit_gy);
    J1 = sum(sum((F-convA(u1)).^2)) + lambda * sum(sum(phi(normgu)));
    u1=u1-step*g1;
    history_g1(i) = norm(g1);
    history_J1(i) = J1;
    history_e1(i) = norm(X-u1);
    
    % denoising
    [gx,gy]=grad_im(u2);
    normgu = sqrt(gx.^2+gy.^2);
    dphigu=dphi(normgu);
    unit_gx = gx./(normgu+1e-7);
    unit_gy = gy./(normgu+1e-7);
    g2 = 2*(u2-F) - lambda*div_champ(dphigu.*unit_gx, dphigu.*unit_gy);
    J2 = sum(sum((F-convA(u2)).^2)) + lambda * sum(sum(phi(normgu)));
    u2=u2-step*g2;
    history_g2(i) = norm(g2);
    history_J2(i) = J2;
    history_e2(i) = norm(X-u2);
end

figure(1);
subplot(311)
plot(1:N, log(history_g1),'DisplayName','deconvolution');
hold on;
plot(1:N, log(history_g2),'DisplayName','denoising');
title('log gradient norm');
legend show;
subplot(312)
plot(1:N, log(history_J1),'DisplayName','deconvolution');
hold on;
plot(1:N, log(history_J2),'DisplayName','denoising');
title('log energy')
legend show;
subplot(313)
plot(1:N, history_e1,'DisplayName','deconvolution');
hold on;
plot(1:N, history_e2,'DisplayName','denoising');
title('error')
legend show;
figure(2);
imshow([X,F,u1,u2],[]);
title('left to right: original, blurred then noise added, Deconvolution, Denoising')
