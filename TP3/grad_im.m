function [gx,gy]=grad_im(im)

%renvoie deux images de même taille que im representant le gradient
%Chcune est une composante du gradient au point correspondant dans l'image 
%suivant les formules du TP.

%taille de l'image (attention la dimension en y est en premier)
[dy,dx]=size(im);

gx=im(:,2:end)-im(:,1:end-1);
%on complète par une colonne de zéros
gx=[gx zeros(dy,1)];

gy=im(2:end,:)-im(1:end-1,:);
gy=[gy;zeros(1,dx)];
