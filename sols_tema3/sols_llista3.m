% Geometria diferencial, FME, UPC
% curs 2018/9
% solucions amb Matlab de problems de la llista 3 
%
% Jaume Amoros, UPC, Barcelona
% 2019/1/16

% variables globals: k=curvatura principal que estiguem seguint (pel
% problema 30)
global k


% problema 12 CURS 2019/20 NO ENTRA
P0=[1 1 1].';
[P,areamin]=fminsearch(@area_estel,P0)
[PH,curvmin]=fminsearch(@curvatura_mitja_vertex,P0)

% problema 30
% k=-0.67; % curvatura principal inicial, calculada a ma o amb linea_curvatura
P0=[1;1;sqrt(5/12)];
% trobem u0,v0 tals que phi(u0,v0)=P per la parametrizacio 'esferica' de
% l'elipsoide
v0=asin(P0(3));
u0=atan2(P0(2)/2,P0(1)/sqrt(3));
w0=[u0;v0];
% busquem les curvatures principals i direccions de curvatura en P
% Es calculen numericament, per diferencies centrades, les derivades 1es i 2es 
% de la parametrizacio de la superficie, perque aquesta superficie es molt
% llisa
% Parametritzacio de la superficie
phi=@(u,v)[sqrt(3)*cos(u).*cos(v);2*sin(u).*cos(v);sin(v)];
% el que ve funciona per qualsevol phi,w0:
% pas de primera derivacio (valor optim en doble precissio)
h=sqrt(1e-15);
% derivades primeres
phiu=(phi(u0+h,v0)-phi(u0-h,v0))/(2*h);
phiv=(phi(u0,v0+h)-phi(u0,v0-h))/(2*h);
% pas de segona derivacio (valor optim en doble precissio) es sqrt(h)
h=sqrt(h);
% derivades segones
phiuu=(phi(u0+h,v0)-2*phi(u0,v0)+phi(u0-h,v0))/(h*h);
phivv=(phi(u0,v0+h)-2*phi(u0,v0)+phi(u0,v0-h))/(h*h);
phiuv=(phi(u0+h,v0+h)-phi(u0-h,v0+h)-phi(u0+h,v0-h)+phi(u0-h,v0-h))/(4*h*h);
% ara ja podem calcular les dues formes fonamentals, el vector normal
% unitari
E=dot(phiu,phiu);
F=dot(phiu,phiv);
G=dot(phiv,phiv);
I=[E,F;F,G];
n=cross(phiu,phiv);
n=n/norm(n);
e=dot(n,phiuu);
f=dot(n,phiuv);
g=dot(n,phivv);
II=[e,f;f,g];
% calculem aplicacio de Weingarten i la diagonalitzem
S=inv(I)*II;
[C,D]=eig(S);
% curvatures mitja i gaussiana en P
H=trace(D)/2
K=det(S)
% curvatures principals
k1=D(1,1)
k2=D(2,2)
% triem seguir la linea de curvatura de k2 perque te el valor absolut mes
% gran
k=k2;
% ara trobem la linea de curvatura integrant el vector unitari de la
% direccio principal de 0 a un final que escollim:
t0=0;
tf=2;
temps=[t0,tf];
% integrem el sistema per la primera linea de curvatura
[t1,uv1]=ode45(@linea_curvatura,temps,w0);
uv1=uv1.';


% i per la segona (k obtinguda com per la primera linea de curvatura)
k=-0.33;
[t2,uv2]=ode45(@linea_curvatura,temps,w0);
uv2=uv2.';


% dibuixem les corbes resultants en la superficie
figure(1)
% la superficie
u=linspace(0,2*pi,61);
v=linspace(-pi/2,pi/2,21);
[u2,v2]=ndgrid(u,v);
x2=sqrt(3)*cos(u2).*cos(v2);
y2=2*sin(u2).*cos(v2);
z2=sin(v2);
mesh(x2,y2,z2);
alpha(0.50)
axis equal

hold on
% linea de curvatura 1
u=uv1(1,:);
v=uv1(2,:);
x=sqrt(3)*cos(u).*cos(v);
y=2*sin(u).*cos(v);
z=sin(v);
plot3(x,y,z,'r','LineWidth',2);
% linea de curvatura 2
u=uv2(1,:);
v=uv2(2,:);
x=sqrt(3)*cos(u).*cos(v);
y=2*sin(u).*cos(v);
z=sin(v);
plot3(x,y,z,'g','LineWidth',2);
hold off

xlabel('x');
ylabel('y');
zlabel('z');
axis equal


    