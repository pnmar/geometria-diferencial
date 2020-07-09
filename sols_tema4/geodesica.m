function jetg=geodesica(t,w)
% Sistema d'edos d'una geodesica en una superficie parametrizada de R^3.
%
% Jaume Amoros, UPC, Barcelona
% 2019/5/3

% en desenvolupament

u=w(1);
v=w(2);
up=w(3);
vp=w(4);

% Problema 15, llista 4: tor que parametrizem com
phi=@(u,v)[(2+cos(u))*cos(v);(2+cos(u))*sin(v);sin(u)];

% Calculem la primera forma fonamental i les seves derivades per derivacio
% numerica (les alternatives mes fines es discuteixen a linea_curvatura.m).
h=sqrt(1e-15);
% derivades primeres
phiu=(phi(u+h,v)-phi(u-h,v))/(2*h);
phiv=(phi(u,v+h)-phi(u,v-h))/(2*h);
% primera forma fonamental
E=dot(phiu,phiu);
F=dot(phiu,phiv);
G=dot(phiv,phiv);
I=[E,F;F,G];

% pas de segona derivacio (valor optim en doble precissio) es sqrt(h)
h=sqrt(h);
% derivades segones
phiuu=(phi(u+h,v)-2*phi(u,v)+phi(u-h,v))/(h*h);
phivv=(phi(u,v+h)-2*phi(u,v)+phi(u,v-h))/(h*h);
phiuv=(phi(u+h,v+h)-phi(u-h,v+h)-phi(u+h,v-h)+phi(u-h,v-h))/(4*h*h);
% derivades de la primera forma fonamental
Eu=2*dot(phiu,phiuu);
Ev=2*dot(phiu,phiuv);
Fu=dot(phiuu,phiv)+dot(phiu,phiuv);
Fv=dot(phiuv,phiv)+dot(phiu,phivv);
Gu=2*dot(phiv,phiuv);
Gv=2*dot(phiv,phivv);

% simbols de Christoffel
Gamma11=I\[Eu/2;Fu-Ev/2];
Gamma12=I\[Ev/2;Gu/2];
Gamma22=I\[Fv-Gu/2;Gv/2];

% sistema d'equacions de la geodesica
jetg=[up;vp;-Gamma11*up*up-2*Gamma12*up*vp-Gamma22*vp*vp];
