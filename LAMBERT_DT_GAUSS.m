function [V1_XYZ, V2_XYZ] = LAMBERT_DT_GAUSS(mu, R1_XYZ, R2_XYZ, dt, isShortWay, tol)
% [V1_XYZ, V2_XYZ] = LAMBERT_DT_GAUSS(mu, R1_XYZ, R2_XYZ, dt, isShortWay)
% Finds the velocities at the start and end of the transfer for
% the arc between R1 and R2 which takes dt seconds. (Vallado 478)
% This method works best for position vectors which are close together
% It will fail for vectors more than 90 degrees apart-
    R1 = norm(R1_XYZ);
    R2 = norm(R2_XYZ);
    
    dTA = atan2(norm(cross(R1_XYZ, R2_XYZ)), dot(R1_XYZ, R2_XYZ));
    if ~isShortWay
        dTA = -dTA;
    end

    l = (R1+R2)/(4*sqrt(R1*R2)*cos(dTA/2)) - 1/2;
    m = (mu * dt^2)/(2*sqrt(R1*R2)*cos(dTA/2))^3;

    y = 1;
    cond = true;
    while cond % do
        x1 = m/y^2 - l;
        x2 = 4/3;
        for i = 0:3
            x2 = x2 + (4/3)*(prod(6:2:6+2*i)*x1^i/prod(5:2:5+2*i));
        end
        y_new = 1 + x2*(l+x1);
        cond = abs(y - y_new) > tol; % while
        y = y_new;
    end

    eccentric = true;
    if eccentric
        dE = 2*acos(1 - 2*x1);
        p = (R1*R2*(1 - cos(dTA)))/(R1 + R2 - 2*sqrt(R1*R2)*cos(dTA/2)*cos(dE/2));
    else
        dH = 2*acosh(1 - 2*x1);
        p = (R1*R2*(1 - cos(dTA)))/(R1 + R2 - 2*sqrt(R1*R2)*cos(dTA/2)*cosh(dH/2));
    end

    [V1_XYZ, V2_XYZ] = KEPLER_FGTA(mu, p, R1_XYZ, R2_XYZ, dTA);
end

