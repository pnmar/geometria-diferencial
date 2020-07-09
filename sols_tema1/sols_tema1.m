% Geometria diferencial
% problemes tema 1

% problema 1
% (b)
t=linspace(-1,1,81);
gammapl=[cosh(t)/4;sinh(t)];
gammam=[-cosh(t)/4;sinh(t)];
plot(gammapl(1,:),gammapl(2,:),'b');
axis equal
hold on
plot(gammam(1,:),gammam(2,:),'r');
% assimptotes
xasimp=[gammam(1,:) gammapl(1,:)];
plot(xasimp,4*xasimp,'k--');
plot(xasimp,-4*xasimp,'k--');
% (c)
theta=linspace(0,2*pi,101);
lambda=(1./(1+tan(theta).^4)).^(1/4);
gamma=[lambda; lambda.*tan(theta)];
figure(2)
plot(gamma(1,:),gamma(2,:));
axis equal
% la solucio es incorrecta quan y<0
% (c) astroide
theta=linspace(0,2*pi,201);
astr=[cos(theta).^3;sin(theta).^3];
figure(3)
plot(astr(1,:),astr(2,:));

% Problema 3
t=linspace(0,3.14,101);
tractr=[sin(t);cos(t)+log(tan(t/2))];
figure(4)
plot(tractr(1,:),tractr(2,:),'LineWidth',2);
axis equal
hold on
plot(tractr(1,1),tractr(2,1),'*g');
plot(tractr(1,end),tractr(2,end),'*r');
hold off

% problema 8
% corba poliedral exemple
t=linspace(0,5,101);
C=[1+t;-t.*t;1+2/3*t.^3];
plot3(C(1,:),C(2,:),C(3,:));
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

longs=arc(C);
anggir=gir(C);

% aproximem la curvatura continua de la corba
pseudocurvatura=anggir./longs(1:end-1);

% busquem centres de curvatura i curvatura
% com a inversa del radi
[kappa,ccurv]=curv_centre(C);

% problema 21
% condicions inicials d'integracio
s0=-1;
s1=1;
P=[0;0;0];
Iden=eye(3);
Y0=[P;Iden(:)];
% integracio del sistema dinamic
[s,Y]=ode45(@funcio_corba,[s0,s1],Y0);

% postproces
Y=Y.';
% corba solucio
gamma=Y(1:3,:);
figure(3)
plot3(gamma(1,:),gamma(2,:),gamma(3,:));
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

% verificacio: es la curvatura 0.5?
longs=arc(gamma);
anggir=gir(gamma);
pseudocurvatura=anggir./longs(1:end-1);

[kappa,ccurv]=curv_centre(gamma);

