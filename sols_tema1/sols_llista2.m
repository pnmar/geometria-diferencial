% Geometria Diferencial, FME, UPC
% 2018/9

% Llista 2: Superficies
% Solucions en Matlab a alguns problemes

% problema 2
% dibuixem l'helicoide
u=linspace(-2*pi,2*pi,41);
v=linspace(-1,1,21);
[u2,v2]=ndgrid(u,v);
x2=sinh(v2).*cos(u2);
y2=sinh(v2).*sin(u2);
z2=u2;
figure(1)
mesh(x2,y2,z2);
axis equal
alpha(0.50) % fraccio de transparencia de la superficie
xlabel('x');
ylabel('y');
zlabel('z');

% problema 3
% dibuix del paraboloide hiperbolic segons paramertizacio (b)
u=linspace(-1.5,1.5,31);
v=linspace(-1,1,21);
[u2,v2]=ndgrid(u,v);
x2=u2+v2;
y2=u2-v2;
z2=4*u2.*v2;
figure(2)
mesh(x2,y2,z2);
% axis equal
alpha(0.50) % fraccio de transparencia de la superficie
xlabel('x');
ylabel('y');
zlabel('z');
% dibuix del paraboloide hiperbolic segons paramertizacio (c)
u=linspace(-1,1,21);
v=linspace(-0.75,0.75,11);
[u2,v2]=ndgrid(u,v);
x2=u2.*cosh(v2);
y2=u2.*sinh(v2);
z2=u2.*u2;
figure(3)
mesh(x2,y2,z2);
% axis equal
alpha(0.50) % fraccio de transparencia de la superficie
xlabel('x');
ylabel('y');
zlabel('z');
