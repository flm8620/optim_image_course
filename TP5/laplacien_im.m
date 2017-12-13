function delta=laplacien_im(im)
%renvoie le laplacien de l'image compatible avec les discretisation 
%utilisées ici

[gx,gy]=grad_im(im);
delta=div_champ(gx,gy);
