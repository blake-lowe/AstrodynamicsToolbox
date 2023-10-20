classdef AnomalyConvert
    %ANOMALYCONVERT Set of static methods for converting between anomalies
    
    methods(Static)
        function [EA] = MA2EA(MA, ECC, tol)
            %[EA] = MA2EA(MA, ECC, tol)
            % Convert Mean Anomaly to Eccentric Anomaly for an elliptical
            % orbit using a Newton-Rhapson solver. (Vallado 65)
            % MA, Mean Anomaly
            % ECC, Eccentricity
            % tol, Numerical tolerance

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

        function [HA] = NA2HA(NA, ECC, tol)
            %[HA] = NA2HA(NA, ECC, tol)
            % Convert Hyperbolic Mean Anomaly to Hyperbolic Eccentric Anomaly 
            % for a hyperbolic orbit using a Newton-Rhapson solver.
            % (Vallado 71)
            % MA, Mean Anomaly
            % ECC, Eccentricity
            % tol, Numerical tolerance

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

        function [BA] = OA2BA(OA)
            %[BA] = OA2BA(OA)
            % Convert Parabolic Mean Anomaly to Parabolic Eccentric Anomaly
            % for a parabolic orbit by solving Barker's Equation (cubic) \
            % using Cardano's Method (Vallado 69, 1027)
            % OA, Parabolic Mean Anomaly
            
            P = 0;
            Q = 3;
            R = -3*OA;
            a = (1/3)*(3*Q - P^2); % = 3
            b = (1/27)*(2*P^3 - 9*P*Q + 27*R); % = R
            del = a^3/27 + b^2/4; % always positive
        
            % Calculate the one real root
            BA = (-b/2 + sqrt(del))^(1/3) + (-b/2 - sqrt(del))^(1/3);
        end

    end
end

