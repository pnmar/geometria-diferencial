% Geometria diferencial, FME, UPC
% curs 2019/20
% solucions amb Matlab de problems de la llista 4 
%
% Jaume Amoros, UPC, Barcelona
% 2020/5/11

% problemes 15,16
% condicio inicial: P=phi(u0,vo)
u0=pi/2;
v0=0;
% condicio inicial w=cos(th0)*phiu+sin(th0)*phiv
% (al agafar-lo amb norma 1 parametrizarem la geodesica per l'arc)
th0=0.55;
up0=cos(th0);
vp0=sin(th0);
jetg0=[u0;v0;up0;vp0];
% interval de temps a integrar el sistema: de t0 a tf
t0=0;
tf=4;
temps=[t0,tf];
% integrem el sistema per la geodesica
[t,jetg]=ode45(@geodesica,temps,jetg0);
jetg=jetg.';

% dibuixem la geodesica en la superficie
figure(1)
% la superficie
nu=61;
nv=81;
u=linspace(0,2*pi,nu);
v=linspace(0,2*pi,nv);
[u2d,v2d]=ndgrid(u,v);
phi2d=[(2+cos(u2d)).*cos(v2d);(2+cos(u2d)).*sin(v2d);sin(u2d)];
x2d=phi2d(1:nu,:);
y2d=phi2d(nu+1:2*nu,:);
z2d=phi2d(2*nu+1:3*nu,:);
mesh(x2d,y2d,z2d);
alpha(0.50)
axis equal
xlabel('x');
ylabel('y');
zlabel('z');

hold on
% geodesica
u=jetg(1,:);
v=jetg(2,:);
% apliquem a u,v la parametritzacio de la superficie per a tenir la
% trajectoria de la geodesica
x=(2+cos(u)).*cos(v);
y=(2+cos(u)).*sin(v);
z=sin(u);
plot3(x,y,z,'r','LineWidth',2);
hold off

xlabel('x');
ylabel('y');
zlabel('z');
axis equal

% calculem la longitud de l'arc de geodesica trobat
% velocitat escalar inicial (que val per tot el domini perque es constant)
h=sqrt(1e-15);
% derivades primeres
phi=@(u,v)[(2+cos(u))*cos(v);(2+cos(u))*sin(v);sin(u)];
phiu0=(phi(u0+h,v0)-phi(u0-h,v0))/(2*h);
phiv0=(phi(u0,v0+h)-phi(u0,v0-h))/(2*h);
% velocitat escalar
velesc=norm(up0*phiu0+vp0*phiv0)
% longitud de l'arc
Lgamma=(temps(2)-temps(1))*velesc



