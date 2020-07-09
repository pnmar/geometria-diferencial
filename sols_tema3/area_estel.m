function area=area_estel(P)
    % vertexs de l'estel, ordenats ciclicament
    V=[1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 1 0 1].';
    % repetim el primer vertex al final de la llista
    Vc=[V,V(:,1)];
    % vectors de P a cada vertex
    Vc=Vc-P*ones(1,size(Vc,2));
    % area
    area=0;
    for k = 1:size(Vc,2)-1
        area = area + norm(cross(Vc(:,k), Vc(:,k+1)));
    end
    area=area/2;
end