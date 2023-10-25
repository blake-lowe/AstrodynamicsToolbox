classdef Frame
    %FRAME Set of static methods for performing coordinate transformations
    
    methods (Static)
        function [out] = body1(in, ANG)
            R = [1 0 0;
                 0 cos(ANG) sin(ANG);
                 0 -sin(ANG) cos(ANG)];
            out = R*in;
        end

        function [out] = body2(in, ANG)
            R = [cos(ANG) 0 -sin(ANG);
                 0 1 0;
                 sin(ANG) 0 cos(ANG)];
            out = R*in;
        end

        function [out] = body3(in, ANG)
            R = [cos(ANG) sin(ANG) 0;
                -sin(ANG) cos(ANG) 0;
                0 0 1];
            out = R*in;
        end

        function [out] = rth2eph(in, TA)
            % Rotation matrix from orbit frame to rotating frame
            R = [cos(TA) sin(TA) 0;
                -sin(TA) cos(TA) 0;
                0 0 1];
            out = R'*in;
        end

        function [out] = rth2xyz(in, AOL, INC, RAAN)
            % in, column vector in orbit-fixed frame to be rotated into inertial frame
            % AOL, argument of latitude (argument of periapsis + true anomaly)
            % INC, inclination
            % RAAN, right ascension of the ascending node
            R = transpose([cos(RAAN)*cos(AOL) - sin(RAAN)*cos(INC)*sin(AOL), sin(RAAN)*cos(AOL) + cos(RAAN)*cos(INC)*sin(AOL), sin(INC)*sin(AOL);
                - cos(RAAN)*sin(AOL) - sin(RAAN)*cos(INC)*cos(AOL), cos(RAAN)*cos(INC)*cos(AOL) - sin(RAAN)*sin(AOL), cos(AOL)*sin(INC);
                sin(RAAN)*sin(INC), -cos(RAAN)*sin(INC), cos(INC)]);
            out = R*in;
        end

        function [out] = rth2bvn(in, FPA)
            out = Frame.body3(in, -FPA);
        end

        function [out] = bvn2rth(in, FPA)
            out = Frame.body3(in, FPA);
        end

        function [out] = eph2rth(in, TA)
            % Rotation matrix from orbit frame to rotating frame
            R = [cos(TA) sin(TA) 0;
                -sin(TA) cos(TA) 0;
                0 0 1];
            out = R*in;
        end

        function [out] = eph2xyz(in, AOP, INC, RAAN)
            % in, column vector in orbit-fixed frame to be rotated into inertial frame
            % AOP, argument of periapsis
            % INC, inclination
            % RAAN, right ascension of the ascending node
            R = transpose([cos(RAAN)*cos(AOP) - sin(RAAN)*cos(INC)*sin(AOP), sin(RAAN)*cos(AOP) + cos(RAAN)*cos(INC)*sin(AOP), sin(INC)*sin(AOP);
                - cos(RAAN)*sin(AOP) - sin(RAAN)*cos(INC)*cos(AOP), cos(RAAN)*cos(INC)*cos(AOP) - sin(RAAN)*sin(AOP), cos(AOP)*sin(INC);
                sin(RAAN)*sin(INC), -cos(RAAN)*sin(INC), cos(INC)]);
            out = R*in;
        end

        function [out] = xyz2eph(in, AOP, INC, RAAN)
            % in, column vector in inertial frame to be rotated into orbit-fixed frame
            % AOP, argument of latitude (argument of periapsis + true anomaly
            % INC, inclination
            % RAAN, right ascension of the ascending node
            R = [cos(RAAN)*cos(AOP) - sin(RAAN)*cos(INC)*sin(AOP), sin(RAAN)*cos(AOP) + cos(RAAN)*cos(INC)*sin(AOP), sin(INC)*sin(AOP);
                - cos(RAAN)*sin(AOP) - sin(RAAN)*cos(INC)*cos(AOP), cos(RAAN)*cos(INC)*cos(AOP) - sin(RAAN)*sin(AOP), cos(AOP)*sin(INC);
                sin(RAAN)*sin(INC), -cos(RAAN)*sin(INC), cos(INC)];
            out = R*in;
        end

        function [out] = xyz2rth(in, AOL, INC, RAAN)
            % in, column vector in inertial frame to be rotated into orbit-fixed frame
            % AOL, argument of latitude (argument of periapsis + true anomaly
            % INC, inclination
            % RAAN, right ascension of the ascending node
            R = [cos(RAAN)*cos(AOL) - sin(RAAN)*cos(INC)*sin(AOL), sin(RAAN)*cos(AOL) + cos(RAAN)*cos(INC)*sin(AOL), sin(INC)*sin(AOL);
                - cos(RAAN)*sin(AOL) - sin(RAAN)*cos(INC)*cos(AOL), cos(RAAN)*cos(INC)*cos(AOL) - sin(RAAN)*sin(AOL), cos(AOL)*sin(INC);
                sin(RAAN)*sin(INC), -cos(RAAN)*sin(INC), cos(INC)];
            out = R*in;
        end


    end
end

