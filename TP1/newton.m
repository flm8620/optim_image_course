function [x_next] = newton(x,G,H)
% x: start point
% d: search step
% J: funtion
% G,H: gradient and hessian
grad = G(x);
hess = H(x);
x_next = x - hess\grad;
end
