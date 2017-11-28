%% test de la conjuguaison des operateurs 
N=256;

ditestx=randn(N,N);
ditesty=randn(N,N);
imtest=randn(N,N);

gx=gradx(imtest);
gy=grady(imtest);

A=sum(gx(:).*ditestx(:))+sum(gy(:).*ditesty(:));
tmp=-div(ditestx,ditesty);
B=sum(tmp(:).*imtest(:));
abs((B-A))/abs(B)

%% creation d'une image de bruit 
imtest=imtest-mean(mean(imtest));
chx0=gradx(imtest);
chy0=grady(imtest);
%% evolution a pas optimal sans preconditionnement 
x0=zeros(N,N);
c=-div(chx0,chy0);
nb=100;
en=zeros(1,nb+1);
for k=1:nb
    en(k)=sum((x0(:)-imtest(:)).^2);
    gx=gradx(x0);
    gy=grady(x0);
    Bx=-div(gx,gy);
    d=-(Bx-c);
    nd=sum(d(:).^2);
    gx=gradx(d);
    gy=grady(d);
    Bd=-div(gx,gy);
    r=sum(sum(d.*Bd));
    t=nd/r;
    x0=x0+t*d;
end
en(end)=sum((x0(:)-imtest(:)).^2);

%% preparation du preconditionnement 
fKx=fft2([1,-1],N,N);
fKy=fft2([1;-1],N,N);
mask=abs(fKx).^2+abs(fKy).^2;
mask(1,1)=1;
maski=1./mask;
maski12=maski.^(0.5);
P= @(y) (real(ifft2(fft2(y).*maski12)));
Pt=@(y)  (real(ifft2(fft2(y).*maski12)));
Pinv=@(y)  (real(ifft2(fft2(y)./maski12)));
%% preconditionnement x=Py on cherche le y optimal
y0=zeros(N,N);
c=-div(chx0,chy0);
c=Pt(c);
nb=100;
enp=zeros(1,nb+1);
for k=1:nb
    x0=P(y0);
    enp(k)=sum((x0(:)-imtest(:)).^2);
    py0=P(y0);
    gx=gradx(py0);
    gy=grady(py0);
    Bx=-div(gx,gy);
    Bx=Pt(Bx);
    d=-(Bx-c);
    nd=sum(d(:).^2);
    Pd=P(d);
    gx=gradx(Pd);
    gy=grady(Pd);
    Bd=-div(gx,gy);
    Bd=Pt(Bd);
    r=sum(sum(d.*Bd));
    t=nd/r;
    y0=y0+t*d;
end
x0=P(y0);
enp(end)=sum((x0(:)-imtest(:)).^2);

