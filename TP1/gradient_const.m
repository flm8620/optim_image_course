function [x_final] = gradient_const(x,G,step)
grad = G(x);
x_final = x-step*grad;
end