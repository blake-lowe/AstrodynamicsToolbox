function [EA] = MA2EA(MA, ECC, tol)
    if (-pi < MA && MA < 0) || (MA > pi)
        EA = MA-ECC;
    else
        EA = MA+ECC;
    end
    
    while 1 
        E_new = EA + (MA-EA+ECC*sin(EA))/(1-ECC*cos(EA));
        if abs(EA - E_new) < tol
            EA = E_new;
            break;
        else
            EA = E_new;
        end
    end
end