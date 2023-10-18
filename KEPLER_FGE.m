function [r1, v1] = KEPLER_FGE(mu, a, dt, dE, r0, v0)
    f = 1 - (a/norm(r0))*(1- cos(dE))
    g = (dt) - sqrt(a^3/mu)*(dE - sin(dE))
    
    r1 = f*r0 + g*v0;
    
    f_dot = -sin(dE)*sqrt(mu*a)/(norm(r0)*norm(r1))
    g_dot = 1 - (a/norm(r1))*(1-cos(dE))
    
    v1 = f_dot*r0 + g_dot*v0;
end

