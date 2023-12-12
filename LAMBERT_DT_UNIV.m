function [V1_XYZ, V2_XYZ] = LAMBERT_DT_UNIV(mu, R1_XYZ, R2_XYZ, dt, isShortWay, tol)
% [V1_XYZ, V2_XYZ] = LAMBERT_DT_UNIV(mu, R1_XYZ, R2_XYZ, dt, isShortWay)
% Finds the velocities at the start and end of the transfer for
% the arc between R1 and R2 which takes dt seconds. (Vallado 492)
% This method will fail if the position vectors have dTA = 180 degrees
% This method may require many iterations
    R1 = norm(R1_XYZ);
    R2 = norm(R2_XYZ);
    
    dTA = atan2(norm(cross(R1_XYZ, R2_XYZ)), dot(R1_XYZ, R2_XYZ));
    A = sqrt(R1*R2*(1 + cos(dTA)));

    if A == 0
        error('LAMBERT_DT_UNIV fails with dTA = 180 degrees.')
    end

    if ~isShortWay
        A = -A;
    end

    psi_n = 0;
    c2 = 1/2;
    c3 = 1/6;
    psi_up = 4*pi^2;
    psi_low = -4*pi;

    cond = true;
    while cond
        y_n = R1 + R2 + A*(psi_n*c3 - 1)/sqrt(c2);
        %if A > 0 && y < 0
            %todo readjust psi_low until y > 0
        %end

        chi_n = sqrt(y_n/c2);
        dt_n = (chi_n^3 * c3 + A * sqrt(y_n))/sqrt(mu);

        if dt_n <= dt
            psi_low = psi_n;
        else
            psi_up = psi_n;
        end

        psi_n1 = (psi_up + psi_low)/2;

        [c2, c3] = PSI2C2C3(psi_n1);

        psi_n = psi_n1;

        cond = abs(real(dt_n) - dt) > tol;
    end

    f = 1 - y_n/R1;
    g_dot = 1 - y_n/R2;
    g = A * sqrt(y_n/mu);
    
    V1_XYZ = (R2_XYZ - f*R1_XYZ)/g;
    V2_XYZ = (g_dot*R2_XYZ - R1_XYZ)/g;
end

