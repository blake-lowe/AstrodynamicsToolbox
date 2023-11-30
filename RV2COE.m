function [SMA, ECC, AOP, INC, RAAN, TA] = RV2COE(MU, R_XYZ, V_XYZ)
    %[SMA, ECC, AOP, INC, RAAN, TA] = RV2COE(MU, R_XYZ, V_XYZ)
    % Converts an inertial position and velocity to Keplerian Orbital Elements
    % (Vallado 113)
    % MU, gravitational parameter in km^2/s^2
    % R_XYZ, position vector in inertial frame
    % V_XYZ, velocity vector in inertial frame

    R = norm(R_XYZ);
    V = norm(V_XYZ);
    H_XYZ = cross(R_XYZ,V_XYZ);
    H = norm(H_XYZ);
    N_XYZ = cross([0,0,1], H_XYZ);
    N = norm(N_XYZ);
    
    % eccentricity vector
    E_XYZ = ((V^2 - MU/R)*R_XYZ - dot(R_XYZ, V_XYZ)*V_XYZ)/MU; 
    ECC = norm(E_XYZ);
    Energy = V^2/2 - MU/R; % specific energy
    
    % semimajor axis and semilatus rectum
    if ECC ~= 1.0
        SMA = -MU/(2*Energy);
        SLR = SMA*(1-ECC^2);
    else
        SMA = inf;
        SLR = H^2/2;
    end
    
    % inclination
    INC = acos(H_XYZ(3)/H);
    
    % right ascension of the ascending node
    RAAN = acos(N_XYZ(1)/N);
    if N_XYZ(2) < 0
        RAAN = 2*pi - RAAN;
    end
    
    % argument of periapsis
    if N == 0 % Elliptical Equatorial
        AOP = acos(E_XYZ(1)/ECC);
        if E_XYZ(2) < 0
            AOP = 2*pi - AOP;
        end
    else
        AOP = acos(dot(N_XYZ, E_XYZ)/(N*ECC));
        if E_XYZ(3) < 0
            AOP = 2*pi - AOP;
        end
    end
    
    % true anomaly
    TA = acos(dot(E_XYZ, R_XYZ)/(ECC*R));
    if dot(R_XYZ, V_XYZ) < 0
        TA = 2*pi - TA;
    end

end