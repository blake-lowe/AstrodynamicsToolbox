function [r1, v1] = KEPLER_FGEA(MU, SMA, dt, dEA, r0, v0)
    f = 1 - (SMA/norm(r0))*(1- cos(dEA))
    g = (dt) - sqrt(SMA^3/MU)*(dEA - sin(dEA))
    
    r1 = f*r0 + g*v0;
    
    f_dot = -sin(dEA)*sqrt(MU*SMA)/(norm(r0)*norm(r1))
    g_dot = 1 - (SMA/norm(r1))*(1-cos(dEA))
    
    v1 = f_dot*r0 + g_dot*v0;
end

