classdef transfer
    %TRANSFER Represents a transfer between two circular orbits
    %For N impulses, Contains 2N Orbits and N Maneuvers
    
    properties
        Orbits
        Maneuvers
        dV
        Time
    end
    
    methods
        function obj = transfer(varargin)
            %TRANSFER Construct an instance of this class
            %transfer('Hohmann', orbit, r_f)
            %transfer('Bi-parabolic', orbit, r_f)
            %transfer('Bi-elliptic', orbit, r_b, r_f)
            %transfer('One-Tangent', orbit, r_f, dTA_b, Periapsis)
            o1 = varargin{2};
            mu = o1.Body.Mu;
            if (o1.ECC ~= 0)
                error('Initial orbit must be circular')
            end

            if strcmp(varargin{1}, 'Hohmann')
                r2 = varargin{3};
                r1 = o1.R;
                a_t = (r1 + r2)/2;

                v1 = sqrt(mu/r1);
                v2 = sqrt(mu/r2);

                vt_a = sqrt((2*mu/r1) - (mu/a_t));
                vt_b = sqrt((2*mu/r2) - (mu/a_t));

                dVa = vt_a - v1;
                dVb = v2 - vt_b;
                m1t = maneuver('alpha', o1, dVa, 0);
                ota = o1.maneuver_execute(m1t);
                otb = ota.propagate_toTA(pi); % Apoapsis
                mt2 = maneuver('alpha', otb, dVb, 0);
                o2 = otb.maneuver_execute(mt2);

                obj.Orbits = [o1, ota, otb, o2];
                obj.Maneuvers = [m1t, mt2];
                obj.dV = dVa + dVb;
                obj.Time = pi*sqrt(a_t^3/mu);
            elseif strcmp(varargin{1}, 'Bi-parabolic')

            elseif strcmp(varargin{1}, 'Bi-elliptic')

            elseif strcmp(varargin{1}, 'One-tangent')
                r2 = varargin{3};
                dTA_b = varargin{4};
                r1 = o1.R;
                R = r1/r2;
                periapsis = varargin{5};
                if periapsis
                    ECC_b = (R-1)/(cos(dTA_b)-R);
                    a_b = r1/(1-ECC_b);
                else
                    %ECC_b = (R-1)/(cos(dTA_b)+R);
                    ECC_b = (1/R-1)/(1/R - cos(dTA_b));
                    a_b = max(r1,r2)/(1+ECC_b);
                end

                v1 = sqrt(mu/r1);
                v2 = sqrt(mu/r2);
                va_b = sqrt((2*mu/r1) - (mu/a_b));
                vb_b = sqrt((2*mu/r2) - (mu/a_b));

                
                
                if periapsis % tangential burn at periapsis
                    dVa = va_b - v1;
                    m1b = maneuver('alpha', o1, dVa, 0);
                    ob_a = o1.maneuver_execute(m1b);
                    ob_b = ob_a.propagate_TA(dTA_b);
                    FPA_b = ob_b.FPA;
                    dVb = sqrt(vb_b^2 + v2^2 - 2*vb_b*v2*cos(FPA_b));
                    alpha_b = pi - asin(sin(FPA_b)*v2/dVb);
                    m2b = maneuver('alpha', ob_b, dVb, -alpha_b);
                    o2 = ob_b.maneuver_execute(m2b);
                    obj.Time = ob_b.TP;
                else % tangential burn at apoapsis
                    FPA_b = atan2(ECC_b*sin(pi-dTA_b), 1+ECC_b*cos(pi-dTA_b)); 
                    dVa = sqrt(va_b^2 + v1^2 - 2*va_b*v1*cos(FPA_b));
                    alpha_b = asin(sin(FPA_b)*va_b/dVa);

                    m1b = maneuver('alpha', o1, dVa, alpha_b);
                    ob_a = o1.maneuver_execute(m1b);
                    ob_b = ob_a.propagate_TA(dTA_b);

                    dVb = v2 - vb_b;

                    m2b = maneuver('alpha', ob_b, dVb, 0);
                    o2 = ob_b.maneuver_execute(m2b);
                    obj.Time = ob_b.TP - ob_a.TP;
                end

                obj.Orbits = [o1, ob_a, ob_b, o2];
                obj.Maneuvers = [m1b, m2b];
                obj.dV = abs(dVa) + abs(dVb);
                    
                %todo
            else
                error('Invalid input mode specified')
            end
        end
    end
end

