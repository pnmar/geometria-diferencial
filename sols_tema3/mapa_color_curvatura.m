% dibuix d'una superficie parametrizada amb mapa de color per a indicar la
% curvatura (Gaussiana o mitja)
%
% Jaume Amoros, UPC, Barcelona
% 2020/5/15

% provat i funciona


% DADES QUE ENTRA L'USUARI
% parametrizacio de la superficie:
% gaussiana de revolucio
sigma2=4;
phi=@(u,v)[u;v;exp((-u.*u-v.*v)/sigma2)]; % vectorialitzar la funcio

% domini de parametrizacio
u0=-5;
uf=5;
v0=-4;
vf=4;

% mallat del domini pel calcul
nu=501;
nv=401;
% FI DE LES DADES QUE ENTRA L'USUARI

% mallem el domini i fem el dibuix pelat de la superficie
% taula de valors en u,v
u1d=linspace(u0,uf,nu);
v1d=linspace(v0,vf,nv);
% taula 2-dimensional, recorre el domini R
[u2d,v2d]=ndgrid(u1d,v1d);
% evaluem vectorialment l'aplicacio de parametrizacio
phi2d=phi(u2d,v2d);
% i ara hem de separar les components x,y,z de phi
x2d=phi2d(1:nu,:);
y2d=phi2d(nu+1:2*nu,:);
z2d=phi2d(2*nu+1:3*nu,:);
% % ja podem dibuixar la superficie (nomes es una prova)
% C2d=cos(u2d.*v2d); % color de prova
% figure(1)
% surf(x2d,y2d,z2d,C2d);
% axis equal
% xlabel('x');
% ylabel('y');
% zlabel('z');
% shading interp
% colorbar
% view(0,90);

% Calcul de les curvatures gaussiana i mitja
% Calcularem les derivades de la parametritzacio numericament, per
% diferencies centrades (hi han dues alternatives si aixo no va be: usar
% derivacio simbolica, o be calcular les derivades 1es i 2es de la
% parametrizacio a ma en cada cas i entrar-les al codi).
% pas de primera derivacio (valor optim en doble precissio)
h=sqrt(1e-15);
% derivades primeres, calculades vectorialment 
% (les columnes de phiu,phiv contenen la corresponent derivada de phi en
% els punts amb un valor concret de v, i cada columna te els vectors 3x1 de
% la derivada de phi per a cada valor concret de u, posats un a sobre de
% l'altre.
phiu=(phi(u2d+h,v2d)-phi(u2d-h,v2d))/(2*h);
phiv=(phi(u2d,v2d+h)-phi(u2d,v2d-h))/(2*h);
% pas de segona derivacio (valor optim en doble precissio) es sqrt(h)
h=sqrt(h);
% derivades segones (en el mateix format que les primeres)
phiuu=(phi(u2d+h,v2d)-2*phi(u2d,v2d)+phi(u2d-h,v2d))/(h*h);
phivv=(phi(u2d,v2d+h)-2*phi(u2d,v2d)+phi(u2d,v2d-h))/(h*h);
phiuv=(phi(u2d+h,v2d+h)-phi(u2d-h,v2d+h)-phi(u2d+h,v2d-h)+phi(u2d-h,v2d-h))/(4*h*h);

% separem les components x,y,z de les derivades (taula rectangular per cada
% una)
xphiu=phiu(1:nu,:);
yphiu=phiu(nu+1:2*nu,:);
zphiu=phiu(2*nu+1:3*nu,:);
xphiv=phiv(1:nu,:);
yphiv=phiv(nu+1:2*nu,:);
zphiv=phiv(2*nu+1:3*nu,:);
% primera forma fonamental (taula rectangular per cada coeficient, per a
% tots els valors de u,v)
E=xphiu.*xphiu+yphiu.*yphiu+zphiu.*zphiu;
F=xphiu.*xphiv+yphiu.*yphiv+zphiu.*zphiv;
G=xphiv.*xphiv+yphiv.*yphiv+zphiv.*zphiv;
% vector normal: com el calculem per tota la taula rectangular de valors de
% u,v alhora, apliquem a pic i pala la formula del producte vectorial per a
% cada component x,y,z del vector normal
xn=yphiu.*zphiv-zphiu.*yphiv;
yn=zphiu.*xphiv-xphiu.*zphiv;
zn=xphiu.*yphiv-yphiu.*xphiv;
% el normalitzem
dA=sqrt(xn.*xn+yn.*yn+zn.*zn);
xn=xn./dA;
yn=yn./dA;
zn=zn./dA;
% segona forma fonamental (novament taula rectangular per a cada
% coeficient)
% primer separem les components x,y,z de les derivades segones
xphiuu=phiuu(1:nu,:);
yphiuu=phiuu(nu+1:2*nu,:);
zphiuu=phiuu(2*nu+1:3*nu,:);
xphiuv=phiuv(1:nu,:);
yphiuv=phiuv(nu+1:2*nu,:);
zphiuv=phiuv(2*nu+1:3*nu,:);
xphivv=phivv(1:nu,:);
yphivv=phivv(nu+1:2*nu,:);
zphivv=phivv(2*nu+1:3*nu,:);
% i ara els coeficients de II
e=xn.*xphiuu+yn.*yphiuu+zn.*zphiuu;
f=xn.*xphiuv+yn.*yphiuv+zn.*zphiuv;
g=xn.*xphivv+yn.*yphivv+zn.*zphivv;
% Calcularem l'aplicacio de Weingarten S=inv(I)*II
% pero novament ho farem vectorialment, calculant cada coeficient de 
% S=[a,b;c,d] per tota la taula rectangular de valors de u,v
detI=dA.*dA; % determinant de I
a=(e.*G-f.*F)./detI;
b=(f.*G-g.*F)./detI;
c=(-e.*F+f.*E)./detI;
d=(-f.*F+g.*E)./detI;
% curvatura Gaussiana (taula rectangular, per tots els valors de u,v)
K=a.*d-b.*c;
% curvatura mitja
H=(a+d)/2;

% dibuixem la superficie colorejant-la segons la curvatura gaussiana
figure(2)
surf(x2d,y2d,z2d,K);
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
shading interp
colorbar
view(0,90);

% dibuixem la superficie colorejant-la segons la curvatura gaussiana
figure(3)
surf(x2d,y2d,z2d,H);
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
shading interp
colorbar
view(0,90);
