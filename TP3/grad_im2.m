function [gx,gy]=grad_im2(im)

%renvoie deux images de m�me taille que im representant le gradient
%Chcune est une composante du gradient au point correspondant dans l'image 
%suivant les formules du TP.

%taille de l'image (attention la dimension en y est en premier)
[dy,dx]=size(im);
gx=im(:,2:end)-im(:,1:end-1);
%on compl�te par une colonne de z�ros
gx=[gx, im(:,end-1)-im(:,end)];

gy=im(2:end,:)-im(1:end-1,:);
gy=[gy;im(end-1,:)-im(end,:)];
