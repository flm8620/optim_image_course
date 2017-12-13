function J=TV(im);
%renvoie la variation totale d'une image
% (phi=t dans le formalisme du TP)

n=norm_grad_im(im);
J=sum(sum(n));
