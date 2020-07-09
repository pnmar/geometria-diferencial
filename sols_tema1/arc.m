function longs = arc(C)
    arestes = C(:, 2:end) - C(:, 1:end-1);
    longs = sqrt(sum(arestes.*arestes));
end