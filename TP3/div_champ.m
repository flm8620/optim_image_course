function d=div_champ(cx,cy)

%cy et cy sont les composantes d'un champ de vecteur
%renvoie la divergence discrète de ce champ de vecteur

[dy,dx]=size(cx);

ddx=zeros(dy,dx);
ddy=ddx;

ddx(:,2:end-1)=cx(:,2:end-1)-cx(:,1:end-2);
ddx(:,1)=cx(:,1);
ddx(:,end)=-cx(:,end-1);


ddy(2:end-1,:)=cy(2:end-1,:)-cy(1:end-2,:);
ddy(1,:)=cy(1,:);
ddy(end,:)=-cy(end-1,:);

d=ddx+ddy;


