function anggir=gir(C)
    % Vector of edges and vector of lengths
    arestes = C(:, 2:end) - C(:, 1:end-1);
    longs = sqrt(sum(arestes.*arestes));
    
    % Angle computation
    prods = sum(arestes(:, 1:end-1).*arestes(:, 2:end));
    cosanggir = prods./(longs(1:end-1).*longs(2:end));
    anggir = acos(cosanggir);
end
