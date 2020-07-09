function v = funcio_corba(s,Y)
    T = Y(4:6);
    N = Y(7:9);
    B = Y(10:12);
    v(1:3) = T;
    v(4:6) = N*curvatura(s);
    v(7:9) = - T*curvatura(s)+B*torsio(s);
    v(10:12) = - N*torsio(s);
    v = v.';
end