function [c2, c3] = PSI2C2C3(psi)
% [c2, c3] = PSI2C2C3(psi)
% For the universal variable formulation find c2 and c3. (Vallado 63)
if psi > 1e-6
    c2 = (1 - cos(sqrt(psi)))/psi;
    c3 = (sqrt(psi) - sin(sqrt(psi)))/sqrt(psi^3);
else
    if psi < -1e-6
        c2 = (1 - cosh(sqrt(-psi)))/psi;
        c3 = (sinh(sqrt(-psi)) - sqrt(-psi))/sqrt((-psi)^3);
    else
        c2 = 1/2;
        c3 = 1/6;
    end
end
end

