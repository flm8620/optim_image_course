function d=div_champ2(cx,cy);

%cy et cy sont les composantes d'un champ de vecteur
%renvoie la divergence discr�te de ce champ de vecteur

[dy,dx]=size(cx);

ddx=zeros(dy,dx);
ddy=ddx;

ddx(:,2:end-1)=cx(:,2:end-1)-cx(:,1:end-2);
ddx(:,1)=zeros(dy,1);
ddx(:,end)=zeros(dy,1);


ddy(2:end-1,:)=cy(2:end-1,:)-cy(1:end-2,:);
ddy(1,:)=zeros(1,dx);
ddy(end,:)=zeros(1,dx);

d=ddx+ddy;


