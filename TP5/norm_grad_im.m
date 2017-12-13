function n=norm_grad_im(im)
%Renvoie la norme du gradient de l'imae en chaque point

[gx,gy]=grad_im(im);
n=(gx.^2+gy.^2).^0.5;
