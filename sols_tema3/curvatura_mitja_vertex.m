function Hp=curvatura_mitja_vertex(P)
% Curvatura mitja en un vertex d'una superficie poliedral a partir del
% seu estel. El vector de curvatura mitja segueix la definicio de M.
% Sullivan, Curvatures of Smooth and Discrete Surfaces.
% Interficie amb l'usuari capada per a facilitar que fminsearch optimitzi
% la funcio.
%
% Jaume Amoros, UPC, Barcelona
% 2019/1/16
%
% provat i funciona


% vertexs de l'estel, ordenats ciclicament
V=[1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 1 0 1].';
% repetim el primer vertex al final de la llista
Vc=[V,V(:,1)];
% vectors normals en les cares de l'estel
Vc=Vc-P*ones(1,size(Vc,2)); % vectors de P a cada vertex
for k=1:size(Vc,2)-1,
    vn=cross(Vc(:,k),Vc(:,k+1));
    vn=vn/norm(vn);
    nu(:,k)=vn;
end;
nu=[nu nu(:,1)];
% vectors de curvatura mitja en cada aresta
for k=1:size(Vc,2)-1,
    H(:,k)=cross(Vc(:,k),nu(:,k)-nu(:,k+1));
end;
% curvatura mitja vector en el vertex 
Hp=sum(H,2)/2;
% curvatura mitja escalar
Hp=norm(Hp);
