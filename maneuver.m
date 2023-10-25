classdef maneuver
    %MANEUVER A representation of an impulsive orbital maneuver
    
    properties
        dV_XYZ
        dV
    end
    
    methods
        function obj = maneuver(varargin)
            %MANEUVER Construct an instance of this class, first term
            %should be a string identifying the input arguments
            %maneuver('XYZ', dV_XYZ)
            %maneuver('BVN', orbit, dV_BVN)
            %maneuver('alpha', orbit, dV, alpha)
            %maneuver('phi_beta, orbit, dV, phi, beta)
            %maneuver('alpha_beta, orbit, dV, alpha, beta)
            if strcmp(varargin{1}, 'XYZ')
                obj.dV_XYZ = varargin{2};
            elseif strcmp(varargin{1}, 'BVN')
                orbit = varargin{2};
                dV_BVN = varargin{3};
                dV_RTH = Frame.bvn2rth(dV_BVN, orbit.FPA);
                obj.dV_XYZ = Frame.rth2xyz(dV_RTH, orbit.AOL, orbit.INC, orbit.RAAN);
            elseif strcmp(varargin{1}, 'alpha')
                orbit = varargin{2};
                dV = varargin{3};
                alpha = varargin{4};
                dV_RTH = [dV*sin(orbit.FPA + alpha); dV*cos(orbit.FPA + alpha); 0];
                obj.dV_XYZ = Frame.rth2xyz(dV_RTH, orbit.AOL, orbit.INC, orbit.RAAN);
            elseif strcmp(varargin{1}, 'phi_beta')
                orbit = varargin{2};
                dV = varargin{3};
                phi = varargin{4};
                beta = varargin{5};
                dV_RTH = dV*[cos(beta)*sin(phi); cos(beta)*cos(phi); sin(beta)];
                obj.dV_XYZ = Frame.rth2xyz(dV_RTH, orbit.AOL, orbit.INC, orbit.RAAN);
            elseif strcmp(varargin{1}, 'alpha_beta')
                orbit = varargin{2};
                dV = varargin{3};
                alpha = varargin{4};
                beta = varargin{5};
                dV_BVN = dV*[cos(beta)*sin(alpha); cos(beta)*cos(alpha); sin(beta)];
                dV_RTH = Frame.bvn2rth(dV_BVN, orbit.FPA);
                obj.dV_XYZ = Frame.rth2xyz(dV_RTH, orbit.AOL, orbit.INC, orbit.RAAN);
            else
                error('Invalid input mode specified')
            end
            obj.dV = norm(obj.dV_XYZ);
        end
        
        function dV_BVN = dV_BVN(obj, orbit)
            %dV_BVN(orbit)
            %Return the maneuver expressed in the passed orbit BVN Frame
            dV_RTH = Frame.xyz2rth(obj.dV_XYZ, orbit.AOL, orbit.INC, orbit.RAAN);
            dV_BVN = Frame.rth2bvn(dV_RTH, orbit.FPA);
        end
    end
end

