function [kappa, ccurv] = curv_centre(C)
    % Edge computation for the given discrete curve
    ar01 = C(:,2:end-1)-C(:,1:end-2);
    ar12 = C(:,3:end)-C(:,2:end-1);
    
    % Scalar products of preceeding and posterior edges
    g11 = sum(ar01.*ar01);
    g12 = sum(ar01.*ar12);
    g22 = sum(ar12.*ar12);
    
    % Curvature center according to given formula
    detG = g11.*g22-g12.*g12;
    lambda1 = (0.5*g11.*g22 - g12.*g12 - 0.5*g12.*g22)./detG;
    lambda2 = (0.5*g11.*g12 + 0.5*g11.*g22)./detG;
    ccurv = C(:,1:end-2) + (ones(3,1)*lambda1).*ar01;
    ccurv = ccurv + (ones(3,1)*lambda2).*ar12;
    
    % curvatura com a invers del radi
    vradi=C(:,1:end-2)-ccurv;
    radi=sqrt(sum(vradi.*vradi));
    kappa=1./radi;
end


