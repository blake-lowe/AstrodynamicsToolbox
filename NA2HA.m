function [HA] = NA2HA(NA, ECC, tol)
    if ECC < 1.6
        if (-pi < NA && NA < 0) || (NA > pi)
            HA = NA-ECC;
        else
            HA = NA+ECC;
        end
    else
        if ECC < 3.6 && abs(NA) > pi
            HA = NA-sign(NA)*ECC;
        else
            HA = NA/(ECC-1);
        end
    end
    
    while 1 %do-while
        H_new = HA + (NA+HA-ECC*sinh(HA))/(ECC*cosh(HA) - 1);

        if abs(HA - H_new) < tol
            HA = H_new;
            break;
        else
            HA = H_new;
        end
    end
end