function [ x ] = minimize( A,b,tau,r, grad, proj )
    zero = zeros(size(A,2),1);
    x = zero; % start from vector 0
    history = zeros(r,1);
    for i = 1:r
        x_old=x;
        x = x-tau*grad(A,b,x);
        x = proj(x);
        history(i)=norm(x-x_old);
    end
    plot(1:r,log(history));
    title('log |x_n-x_{n-1}|');
end

