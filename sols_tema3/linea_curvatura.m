function dirpral=linea_curvatura(t,w)
global k % la curvatura principal que estem seguint
% punt en que volem la direccio de curvatura principal
u=w(1);
v=w(2);
% parametrizacio de la superficie (s'ha d'entrar a ma per cada superficie
% on es faci el calcul)
% Problema 30, llista 3: elipsoide que parametrizem com
phi=@(u,v)[sqrt(3)*cos(u)*cos(v);2*sin(u)*cos(v);sin(v)];

% Es calculen numericament, per diferencies centrades, les derivades 1es i 2es 
% de la parametrizacio de la superficie. Aixo donara error apreciable si
% les segones derivades de phi oscil.len molt.
% Alternatives: 
% 1. Que Matlab les derivi simbolicament en el main i les guardi en
% variables globals (programacio mes complicada, incompatible amb Octave).
% 2. Entra aqui a ma el valor analitic de cada derivada primera o segona (solucio
% senzilla i mes acurada, pero porta molta mes feina plantejar cada calcul)
% pas de primera derivacio (valor optim en doble precissio)
h=sqrt(1e-15);
% derivades primeres
phiu=(phi(u+h,v)-phi(u-h,v))/(2*h);
phiv=(phi(u,v+h)-phi(u,v-h))/(2*h);
% pas de segona derivacio (valor optim en doble precissio) es sqrt(h)
h=sqrt(h);
% derivades segones
phiuu=(phi(u+h,v)-2*phi(u,v)+phi(u-h,v))/(h*h);
phivv=(phi(u,v+h)-2*phi(u,v)+phi(u,v-h))/(h*h);
phiuv=(phi(u+h,v+h)-phi(u-h,v+h)-phi(u+h,v-h)+phi(u-h,v-h))/(4*h*h);

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

% controlem quina linea de curvatura estavem seguint: es la que te
% curvatura principal propera a la que teniem
vaps=diag(D); % descomentar aquesta linea per a veure les dues curvatures principals en el punt
[var,bo]=min(abs(vaps-k)); 
dirpral=C(:,bo); % el vector propi ja ve normalitzat, direccio que hem de seguir
k=vaps(bo); % actualitzem la curvatura que seguim

