function [x_final] = gradient_optim(x,G,H)
% x: start point
% d: search step
% J: funtion
% G,H: gradient and hessian
grad = G(x);
[x_final,~] = newton_search(x,-grad,G,H);
end

