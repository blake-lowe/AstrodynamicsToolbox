function [E] = MA2EA(M, e, tol)
    if (-pi < M && M < 0) || (M > pi)
        E = M-e;
    else
        E = M+e;
    end
    
    while 1 
        E_new = E + (M-E+e*sin(E))/(1-e*cos(E));
        if abs(E - E_new) < tol
            E = E_new;
            break;
        else
            E = E_new;
        end
    end
end