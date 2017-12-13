function [u,ener]=Chambolle(im,mu,nb,tau)

[dy,dx]=size(im);
px=zeros(dy,dx);
py=px;
immu=im/mu;

ener(1)=TV(im);

for k=1:nb
    [gx,gy]=grad_im(div_champ(px,py)-immu);
    numer=1+tau*((gx.^2+gy.^2).^0.5);
    px=(px+tau*gx)./numer;
    py=(py+tau*gy)./numer;
    u=im-mu*div_champ(px,py);
    %ener(k+1)=TV(u)+1/(2*mu)*sum((im(:)-u(:)).^2);
    
end
