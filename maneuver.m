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

        function dV_RTH = dV_RTH(obj, orbit)
            %dV_RTH(orbit)
            %Return the maneuver expressed in the passed orbit RTH Frame
            dV_RTH = Frame.xyz2rth(obj.dV_XYZ, orbit.AOL, orbit.INC, orbit.RAAN);
        end

        function dV_EPH = dV_EPH(obj, orbit)
            %dV_EPH(orbit)
            %Return the maneuver expressed in the passed orbit EPH Frame
            dV_EPH = Frame.xyz2eph(obj.dV_XYZ, orbit.AOP, orbit.INC, orbit.RAAN);
        end

        function dV_BVN = dV_BVN(obj, orbit)
            %dV_BVN(orbit)
            %Return the maneuver expressed in the passed orbit BVN Frame
            dV_RTH = Frame.xyz2rth(obj.dV_XYZ, orbit.AOL, orbit.INC, orbit.RAAN);
            dV_BVN = Frame.rth2bvn(dV_RTH, orbit.FPA);
        end

        function obj = reverse(obj)
            obj.dV_XYZ = -obj.dV_XYZ;
        end

        function obj = combine(obj, maneuver2)
            obj.dV_XYZ = obj.dV_XYZ + maneuver2.dV_XYZ;
            obj.dV = norm(obj.dV_XYZ);
        end

        function plot_EPH_tikz(obj, orbit)
            vs = 2;
            rs = 1e4;
            r = orbit.R_EPH/rs;
            v1 = orbit.V_EPH/vs;
            dv = obj.dV_EPH(orbit)/vs;
            v2 = v1 + dv;
            gamma1 = orbit.FPA;
            gamma2 = orbit.maneuver_execute(obj).FPA;
            alpha = atan2(norm(cross(dv, v1)), dot(dv, v1));
            alpha = rad2deg(alpha);

            fprintf('\\begin{figure}[htpb]\n')
            fprintf('\\centering\n')
            fprintf('\\fbox{\n')
            fprintf('\\begin{tikzpicture}[scale=1]\n')
            fprintf('\\coordinate (sc) at (%0.4f, %0.4f); \n', r(1), r(2))
            fprintf('\\node at (0,0) {$\\%s$};\n', orbit.Body.Name)
            fprintf('\\draw [thin, dotted] (0,0) -- (sc);\n')
            fprintf('\\draw [->] (0,0) -- (0,1) node [above] {$\\hat p$};\n')
            fprintf('\\draw [->] (0,0) -- (1,0) node [above] {$\\hat e$};\n')
            e_hat = [1;0;0];
            p_hat = [0;1;0];
            r_hat = Frame.rth2eph(e_hat, orbit.TA);
            t_hat = Frame.rth2eph(p_hat, orbit.TA);
            fprintf('\\draw [thin, ->] (sc) -- ($(sc) + (%0.4f,%0.4f)$) node [left] {$\\hat r$};\n', r_hat(1), r_hat(2))
            fprintf('\\draw [->] (sc) -- ($(sc) + (%0.4f,%0.4f)$) node [below right] {$\\hat \\theta$};\n', t_hat(1), t_hat(2))
            fprintf('\\draw [dashed] (sc) -- ($(sc) - (%0.4f,%0.4f)$) node [right] {l.h.};\n', t_hat(1), t_hat(2))
            fprintf('\\draw [fill=black] (sc) circle [radius = 0.025];\n')
            fprintf('\\draw [thick, ->, blue] (sc) -- ($(sc) + (%0.4f, %0.4f)$) node [above right] {$\\bar v^-$};\n', v1(1), v1(2))
            fprintf('\\draw [dashed, blue] (sc) -- ($(sc) + 1.2*(%0.4f, %0.4f)$);\n', v1(1), v1(2))
            fprintf('\\draw [thick, ->, red] (sc) -- ($(sc) + (%0.4f, %0.4f)$) node [below right] {$\\bar v^+$};\n', v2(1), v2(2))
            fprintf('\\draw [thick, ->, orange] ($(sc) + (%0.4f, %0.4f)$) -- ($(sc) + (%0.4f, %0.4f)$) node [below left] {$\\Delta \\bar v$};\n', v1(1), v1(2), v2(1), v2(2))
            a0 = rad2deg(orbit.TA);
            at = rad2deg(orbit.TA+pi/2);
            a1 = rad2deg(orbit.TA+pi/2-gamma1);
            a2 = rad2deg(orbit.TA+pi/2-gamma2);
            fprintf('\\draw[thick,magenta, |->] ([shift=(0:0.25)]0,0) arc (0:%0.4f:0.25) node [midway, above right] {$\\theta^*$};\n', a0)
            fprintf('\\draw[thick,magenta, |->] ([shift=(%0.4f:0.25)]sc) arc (%0.4f:%0.4f:0.25) node [above] {$\\gamma^-$};\n', at, at, a1)
            fprintf('\\draw[thick,magenta, |->] ([shift=(%0.4f:0.75)]sc) arc (%0.4f:%0.4f:0.75) node [midway, left] {$\\gamma^+$};\n', at, at, a2)
            fprintf('\\draw[thick,magenta, |->] ([shift=(%0.4f:1.25)]sc) arc (%0.4f:%0.4f:1.25) node [midway, left] {$\\Delta \\gamma$};\n',  a1, a1, a2)
            fprintf('\\draw[thick,magenta, |->] ([shift=(%0.4f:0.25)]$(sc) + (%0.4f, %0.4f)$) arc (%0.4f:%0.4f:0.25) node [midway, below left] {$\\alpha$};\n', a1, v1(1), v1(2), a1, a1 - alpha)
            fprintf('\\draw[thick,cyan, <->] ([shift=(%0.4f:0.25)]$(sc) + (%0.4f, %0.4f)$) arc (%0.4f:%0.4f:0.25) node [midway, below right] {$\\eta$};\n', at+180, v1(1), v1(2), at+180, a1-alpha)
            fprintf('\\end{tikzpicture}\n')
            fprintf('}\n')
            fprintf('\\end{figure}\n')
        end
    end
end

