function [v0_XYZ, v_XYZ] = KEPLER_FGTA(mu, p, r0_XYZ, r_XYZ, dTA)
    r0 = norm(r0_XYZ);
    r = norm(r_XYZ);
    f = 1 - (r/p)*(1 - cos(dTA));
    g = (r*r0*sin(dTA))/(sqrt(mu)*p);
    f_dot = sqrt(1/p)*tan(dTA)*((1 - cos(dTA))/p - 1/r - 1/r0);
    g_dot = 1 - (r0/p)*(1 - cos(dTA));
    
    v0_XYZ = (r_XYZ - f*r0_XYZ)/g;
    %v_XYZ = (g_dot*r_XYZ - r0_XYZ)/g;
    v_XYZ = (f_dot*r0_XYZ + g_dot*v0_XYZ);
end

