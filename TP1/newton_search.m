function [x_final,t] = newton_search(x,d,G,H)
% x: start point
% d: search step
% J: funtion
% G,H: gradient and hessian
t=0;
for i=1:50
    grad = G(x+d*t);
    hess = H(x+d*t);
    g = grad'*d;
    h = d'*hess*d;
    t = t - g/h;
    %fprintf('t=%f\n',t);
end
x_final = x+d*t;
end

