function [H] = NA2HA(N, e, tol)
    if e < 1.6
        if (-pi < N && N < 0) || (N > pi)
            H = N-e;
        else
            H = N+e;
        end
    else
        if e < 3.6 && abs(N) > pi
            H = N-sign(N)*e;
        else
            H = N/(e-1);
        end
    end
    
    while 1 %do-while
        H_new = H + (N+H-e*sinh(H))/(e*cosh(H) - 1);

        if abs(H - H_new) < tol
            H = H_new;
            break;
        else
            H = H_new;
        end
    end
end