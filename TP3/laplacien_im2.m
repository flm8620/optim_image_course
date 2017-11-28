function delta=laplacien_im2(im)
%renvoie le laplacien de l'image compatible avec les discretisation 
%utilis�es ici

[gx,gy]=grad_im2(im);
delta=div_champ2(gx,gy);
