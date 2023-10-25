classdef transfer
    %TRANSFER Represents a transfer between two circular orbits
    %Contains N Orbits and N-1 Maneuvers
    
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
            %transfer('Bi-elliptic', orbit, r_b, r_f)
            %transfer('One-Tangent', orbit, r_f, TA_b)
            o1 = varargin{2};
            if (o1.ECC ~= 0)
                error('Initial orbit must be circular')
            end

            if strcmp(varargin{1}, 'Hohmann')
                r2 = varargin{3};
                r1 = o1.R;
                a_t = (r1 + r2)/2;

                mu = o1.Body.Mu;
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

                obj.Orbits = [o1, otb, o2];
                obj.Maneuvers = [m1t, mt2];
                obj.dV = dVa + dVb;
                obj.Time = pi*(a_t^3/mu);

            elseif strcmp(varargin{1}, 'Bi-elliptic')

            elseif strcmp(varargin{1}, 'One-Tangent')

            else
                error('Invalid input mode specified')
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

