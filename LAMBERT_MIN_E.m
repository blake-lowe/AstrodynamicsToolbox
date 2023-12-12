function [a_min, e_min, t_min, V1_XYZ] = LAMBERT_MIN_E(mu, R1_XYZ, R2_XYZ, isShortWay)
% [a_min, e_min, t_min, V1_XYZ] = LAMBERT_MIN_E(R1_XYZ, R2_XYZ, isShortWay)
% Finds the orbital elements for the minimum energy transfer
% arc between R1 and R2. (Vallado 475)

    R1 = norm(R1_XYZ);
    R2 = norm(R2_XYZ);
    
    dTA = atan2(norm(cross(R1_XYZ, R2_XYZ)), dot(R1_XYZ, R2_XYZ));
    c = sqrt(R1^2 + R2^2 - 2*R1*R2*cos(dTA));
    s = (R1 + R2 + c)/2;
    a_min = s/2;
    p_min = (R1*R2/c)*(1-cos(dTA));
    e_min = sqrt(1 - 2*p_min/s);
    alpha = pi;
    beta = 2*asin(sqrt((s-c)/s));
    if isShortWay
        t_min = sqrt(a_min^3/mu)*(alpha - (beta - sin(beta)));
    else
        t_min = sqrt(a_min^3/mu)*(alpha + (beta - sin(beta)));
    end
    
    % use f and g functions to get velocity
    V1_XYZ = (sqrt(mu*p_min)/(R1*R2*sin(dTA)))*(R2_XYZ - (1 - (R2/p_min)*(1 - cos(dTA)))*R1_XYZ);
    
    % reverse velocity if long way
    if ~isShortWay
        V1_XYZ = -V1_XYZ;
    end

end

