classdef orbit
    %ORBIT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties % Distances in km, Angles in rad, Time in s
        Body        % Primary body, as a Capitalized string
        MU          % Gravitational parameter of Body
        SMA         % Semi-major axis
        ECC         % Eccentricity
        AOP         % Argument of Periapsis
        INC         % Inclination
        RAAN        % Right Ascencion of the Ascending Node
        TA          % True Anomaly
        AOL         % Argument of Latitude
        Ascending   % Ascending or Descending, as a logical
        Type        % Ellipse(0), Parabola(1), Hyperbola(2), as an int
        RP          % Radius at Periapsis
        RA          % Radius at Apoapsis
        EA          % Eccentric Anomaly
        MA          % Mean Anomaly
        BA          % Parabolic Anomaly
        OA          % Parabolic Mean Anomaly
        HA          % Hyperbolic Anomaly
        NA          % Hyperbolic Mean Anomaly
        MM          % Mean Motion
        XA          % Universal Mean Anomaly
        Period      % Period
        TP          % Time since periapsis
        SLR           % Semi-latus Rectum
        Energy      % Specific Energy
        R           % Magnitude of Position Vector
        V           % Magnitude of Velocity Vector
        H           % Magnitude of Angular Momentum Vector
        FPA         % Flight Path Angle relative to local horizon
        R_RTH       % Position vector in rotating orbit frame
        V_RTH       % Velocity vector in rotating orbit frame
        R_EPH       % Position vector in fixed orbit frame
        V_EPH       % Velcotiy vector in fixed orbit frame
        R_XYZ       % Position vector in inertial frame
        V_XYZ       % Velocity vector in inertial frame
        H_XYZ       % Angular Momentum Vector in inertial frame
        E_XYZ       % Eccentricity Vector in inertial frame
        N_XYZ       % Nodal Vector in inertial frame
    end
    
    methods
        function obj = orbit(varargin)
            %ORBIT Construct an instance of this class, first term should
            %be a string indentifying the input arguments
            if strcmp(varargin{1},'COE') %Keplerian Orbital Elements
                obj.Body = varargin{2};
                obj.SMA = varargin{3};
                obj.ECC = varargin{4};
                obj.AOP = varargin{5};
                obj.INC = varargin{6};
                obj.RAAN = varargin{7};
                obj.TA = varargin{8};
            end

            if strcmp(varargin{1},'RARP')
                obj.Body = varargin{2};
                RA = varargin{3};
                RP = varargin{4};
                obj.SMA = (RA + RP)/2;
                obj.ECC = (RA - RP)/(RA + RP);
                obj.AOP = varargin{5};
                obj.INC = varargin{6};
                obj.RAAN = varargin{7};
                obj.TA = varargin{8};
            end

            if strcmp(varargin{1}, 'RV') %Position and velocity
                obj.Body = varargin{2};
                [obj.SMA, obj.ECC, obj.AOP, obj.INC, obj.RAAN, obj.TA] =...
                    RV2COE(Body_Grav_param(obj.Body), varargin{3}, varargin{4});
            end

            if strcmp(varargin{1}, 'RRP') %Two Positions and P
                obj.Body = varargin{2};
                % TODO using f and g functions to get v1 then RV2COE(r1,v1)

            end

            obj = obj.update_properties();
        end

        function obj = propagate_TA(obj, dTA)
            obj.TA = obj.TA + dTA;
            obj = obj.update_properties();
        end

        function obj = propagate_toTA(obj, TA)
            if TA > 2*pi || (TA < -pi)
                error('Invalid New TA')
            end
            if obj.TA > TA && obj.ECC < 1
                TA = TA + 2*pi;
            else
                error('New TA must be greater than current TA for non-elliptical orbits')
            end

            dTA = TA - obj.TA;
            obj = obj.propagate_TA(dTA);
        end

        function obj = propagate_Time(obj, dt)
            obj.XA = obj.XA + dt*obj.MM;
            if obj.Type == 0 % Ellipse
                obj.MA = obj.XA;
                obj.EA = MA2EA(obj.MA, obj.ECC, 1e-12);
                obj.TA = mod(2*atan(sqrt((1+obj.ECC)/(1-obj.ECC))*tan(obj.EA/2)), 2*pi);
            elseif obj.Type == 1 % Parabola
                obj.OA = obj.XA;
                obj.BA = OA2BA(obj.OA);
                obj.TA = 2*atan(obj.BA);
            else % Hyperbola
                obj.NA = obj.XA;
                obj.HA = NA2HA(obj.NA, obj.ECC, 1e-12);
                obj.TA = mod(2*atan(sqrt((1+obj.ECC)/(1-obj.ECC))*tanh(obj.HA/2)), 2*pi);
            end

            obj = obj.update_properties();
        end

        function [sma, ecc, aop, inc, raan, ta] = maneuver(obj, dV, alpha)
            % TODO maneuver

            sma = obj.SMA;
            ecc = obj.ECC;
            aop = obj.AOP;
            inc = obj.INC;
            raan = obj.RAAN;
            ta = obj.TA;
        end

        function [sma, ecc, aop, inc, raan, ta] = maneuver_XYZ(obj, dV)
            if length(dV) ~= 3
                error('dV must be have length 3')
            end
            % TODO maneuver

            sma = obj.SMA;
            ecc = obj.ECC;
            aop = obj.AOP;
            inc = obj.INC;
            raan = obj.RAAN;
            ta = obj.TA;
        end

        function obj = update_properties(obj)
            %UPDATE_PROPERTIES, calculates all derived non COE properties

            %MU          % Gravitational parameter of Body
            obj.MU = Body_Grav_Param(obj.Body);

            % P           % Semi-latus Rectum
            obj.SLR = obj.SMA*(1 - obj.ECC^2);

            %AOL         % Argument of Latitude
            obj.AOL = obj.AOP + obj.TA;

            % Ascending   % Ascending or Descending, as a logical
            obj.Ascending = (-pi < obj.TA && obj.TA < 0) || (obj.TA > pi);

            % Type        % Ellipse, Parabola, Hyperbola, as a string
            % EA,MA,BA,OA,HA,NA        % * Anomalies
            % MM          % Mean Motion
            % RP          % Radius at Periapsis
            % RA          % Radius at Apoapsis
            obj.EA = NaN;
            obj.MA = NaN;
            obj.BA = NaN;
            obj.OA = NaN;
            obj.HA = NaN;
            obj.NA = NaN;
            if obj.ECC < 1 % Ellipse
                obj.Type = 0;
                % EA          % Eccentric Anomaly
                obj.EA = mod(2*atan(sqrt((1-obj.ECC)/(1+obj.ECC))*tan(obj.TA/2)), 2*pi);
                % MA          % Mean Anomaly
                obj.MA = obj.EA - obj.ECC*sin(obj.EA);
                obj.MM = sqrt(obj.MU/obj.SMA^3);
                obj.XA = obj.MA;
                obj.RP = obj.SMA*(1 - obj.ECC);
                obj.RA = obj.SMA*(1 + obj.ECC);
            elseif obj.ECC == 1 % Parabola
                obj.Type = 1;
                if obj.TA > pi
                    obj.TA = obj.TA - 2*pi;
                end
                % BA          % Parabolic Anomaly
                obj.BA = tan(obj.TA/2);
                % OA          % Parabolic Mean Anomaly
                obj.OA = obj.BA + obj.BA^3/3;
                obj.MM = 2*sqrt(obj.MU/obj.SLR^3);
                obj.XA = obj.OA;
                obj.RP = obj.SLR/2;
                obj.RA = Inf;
            else % Hyperbola
                obj.Type = 2;
                if obj.TA > pi
                    obj.TA = obj.TA - 2*pi;
                end
                % HA          % Hyperbolic Anomaly
                obj.HA = mod(2*atanh(sqrt((1-obj.ECC)/(1+obj.ECC))*tan(obj.TA/2)), 2*pi);
                % NA          % Hyperbolic Mean Anomaly
                obj.NA = obj.ECC*sinh(obj.HA) - obj.HA;
                obj.MM = sqrt(obj.MU/(-obj.SMA^3));
                obj.XA = obj.NA;
                obj.RP = obj.SMA*(1-obj.ECC);
                obj.RA = NaN;
            end
            
            % Period      % Period
            obj.Period = 2*pi/obj.MM;

            % TP          % Time since periapsis
            obj.TP = obj.XA/obj.MM;

            % Energy      % Specific Energy
            if obj.Type == 1 % parabola
                obj.Energy = 0;
            else
                obj.Energy = -obj.MU/(2*obj.SMA);
            end

            % R       % Magnitude of Position Vector
            obj.R = obj.SLR/(1 + obj.ECC*cos(obj.TA));

            % V       % Magnitude of Velocity Vector
            obj.V = sqrt(obj.MU*(2/obj.R - 1/obj.SMA));

            % H       % Magnitude of Angular Momentum Vector
            obj.H = sqrt(obj.MU*obj.SLR);

            % FPA         % Flight Path Angle relative to local horizon
            obj.FPA = acos(obj.H/(obj.R*obj.V));
            if obj.Ascending
                obj.FPA = -obj.FPA;
            end

            % R_RTH       % Position vector in rotating orbit frame
            obj.R_RTH = [obj.R; 0; 0];
            % V_RTH       % Velocity vector in rotating orbit frame
            obj.V_RTH = [obj.V*sin(obj.FPA); obj.V*cos(obj.FPA); 0];

            % R_EPH       % Position vector in fixed orbit frame
            obj.R_EPH = rth2eph(obj.R_RTH, obj.TA);
            % V_EPH       % Velcotiy vector in fixed orbit frame
            obj.V_EPH = rth2eph(obj.V_RTH, obj.TA);

            % R_XYZ       % Position vector in inertial frame
            obj.R_XYZ = rth2xyz(obj.R_RTH, obj.AOL, obj.INC, obj.RAAN);
            % V_XYZ       % Velocity vector in inertial frame
            obj.V_XYZ = rth2xyz(obj.V_RTH, obj.AOL, obj.INC, obj.RAAN);

            % H_XYZ       % Angular Momentum Vector in inertial frame
            obj.H_XYZ = cross(obj.R_XYZ, obj.V_XYZ);

            % E_XYZ       % Eccentricity Vector in inertial frame
            obj.E_XYZ = (1/obj.MU)*((obj.V^2 - obj.MU/obj.R)*obj.R_XYZ - dot(obj.R_XYZ, obj.V_XYZ)*obj.V_XYZ);

            % N_XYZ       % Nodal Vector in inertial frame
            obj.N_XYZ = cross([0;0;1], obj.H_XYZ);
            obj.N_XYZ = obj.N_XYZ./norm(obj.N_XYZ);
        end
        
        function get_degrees(obj)
            % TODO print all angle measurements in degrees

            % TODO find out desired units for HA
        end

        function plot_EPH(obj, n, draw_body, draw_pos, draw_vecs, draw_apsides, draw_nodes, draw_aux, legend_labels) % Only works for ellipses for now
            ta_vec = linspace(0, 2*pi, n);
            r_mag_vec = obj.SLR./(1+obj.ECC*cos(ta_vec));
            r_e_vec = r_mag_vec.*cos(ta_vec);
            r_p_vec = r_mag_vec.*sin(ta_vec);
            aux_e_vec = -obj.SMA*obj.ECC + obj.SMA*cos(ta_vec);
            aux_p_vec = obj.SMA*sin(ta_vec);
            
            figure()
            hold on
            plot(r_e_vec, r_p_vec)
            xlabel("$$\hat{e}$$ direction [km]",'Interpreter','Latex')
            ylabel("$$\hat{p}$$ direction [km]",'Interpreter','Latex')
            title("Orbit Plot, $$\hat{e},\hat{p},\hat{h}$$ Frame, AAE 532 Blake Lowe", 'Interpreter','Latex')
            axis equal
            if (draw_aux)
                plot(aux_e_vec, aux_p_vec)
            end
            %scatter(0, 0, 'black')
            scatter(0, 0, 'black', 'x') % focus
            scatter(-obj.SMA*obj.ECC, 0, 'black', '+') %center
            if (draw_body)
                p = nsidedpoly(1000, 'Center', [0 0], 'Radius', Body_Radius(obj.Body));
                plot(p, 'FaceColor', 'b')
            end
            if (draw_apsides)
                plot([-obj.RA,obj.RP], [0,0], 'black--')
            end
            if (draw_nodes)
                % TODO
                %N_EPH = xyz2eph(N_XYZ, obj.AOP, obj.INC, obj.RAAN);
                plot([-cos(-obj.AOP), cos(-obj.AOP)]*obj.SMA, [-sin(-obj.AOP), sin(-obj.AOP)]*3*obj.SMA, 'red--')
            end
            if (draw_pos)
                scatter(obj.R_EPH(1), obj.R_EPH(2), 'black')
            end
            if (draw_vecs)
                % TODO unit vectors, vectors, and angles
                plot([0,obj.R_EPH(1)], [0,obj.R_EPH(2)], 'black')
            end
            ylim([min(r_p_vec) - 0.1*obj.SMA, max(r_p_vec) + 0.1*obj.SMA])

            if exist('legend_labels', 'var')
                legend(legend_labels)
            end

            hold off
        end

        function plot_EPH_overlay(obj, n, orbit2, draw_body, draw_pos, draw_vecs, draw_apsides, legend_labels)
            dAOP = orbit2.AOP - obj.AOP;

            ta_vec = linspace(0, 2*pi, n);
            r_mag_vec = obj.SLR./(1+obj.ECC*cos(ta_vec));
            r_e_vec = r_mag_vec.*cos(ta_vec);
            r_p_vec = r_mag_vec.*sin(ta_vec);

            r2_mag_vec = orbit2.SLR./(1+orbit2.ECC*cos(ta_vec + dAOP));
            r2_e_vec = r2_mag_vec.*cos(ta_vec);
            r2_p_vec = r2_mag_vec.*sin(ta_vec);

            figure()
            hold on
            plot(r_e_vec, r_p_vec)
            plot(r2_e_vec, r2_p_vec)
            xlabel("$$\hat{e}$$ direction [km]",'Interpreter','Latex')
            ylabel("$$\hat{p}$$ direction [km]",'Interpreter','Latex')
            title("Overlaid Orbit Plot, $$\hat{e},\hat{p},\hat{h}$$ Frame, AAE 532 Blake Lowe", 'Interpreter','Latex')
            axis equal

            %scatter(0, 0, 'black')
            scatter(0, 0, 'black', 'x') % focus
            scatter(-obj.SMA*obj.ECC, 0, 'black', '+') %center
            if (draw_body)
                p = nsidedpoly(1000, 'Center', [0 0], 'Radius', Body_Radius(obj.Body));
                plot(p, 'FaceColor', 'b')
            end
            if (draw_apsides)
                plot([-obj.RA,obj.RP], [0,0], 'black--')
                plot([-orbit2.RA,orbit2.RP]*cos(-dAOP), [-orbit2.RA,orbit2.RP]*sin(-dAOP), 'black--')
            end
            if (draw_pos)
                scatter(obj.R_EPH(1), obj.R_EPH(2), 'black')
            end
            if (draw_vecs)
                % TODO unit vectors, vectors, and angles
                plot([0,obj.R_EPH(1)], [0,obj.R_EPH(2)], 'black')
            end
            if orbit2.SMA > obj.SMA
                ylim([-1.1, 1.1]*orbit2.SMA*(1 + sin(dAOP)/2))
            else
                ylim([-1.1, 1.1]*obj.SMA*(1 + sin(dAOP)/2))
            end

            ylim([min(min(r_p_vec), min(r2_p_vec)) - 0.1*max(obj.SMA, orbit2.SMA), ...
                  max(max(r_p_vec), max(r2_p_vec)) + 0.1*max(obj.SMA, orbit2.SMA)])

            if exist('legend_labels', 'var')
                legend(legend_labels)
            end

            hold off
        end

        function plot_XYZ(obj, n, draw_body, draw_pos)
            ta_vec = linspace(0, 2*pi, n);
            r_mag_vec = obj.SLR./(1+obj.ECC*cos(ta_vec));
            r_e_vec = r_mag_vec.*cos(ta_vec);
            r_p_vec = r_mag_vec.*sin(ta_vec);
            r_eph_vec = [r_e_vec; r_p_vec; zeros(1,n)];
            r_xyz_vec = zeros(3, n);
            for i = 1:n
                r_xyz_vec(:,i) = eph2xyz(r_eph_vec(:,i), obj.AOP, obj.INC, obj.RAAN);
            end
            r_x_vec = r_xyz_vec(1,:);
            r_y_vec = r_xyz_vec(2,:);
            r_z_vec = r_xyz_vec(3,:);

            figure()
            hold on
            plot3(r_x_vec, r_y_vec, r_z_vec)
            title("Orbit Plot, $$\hat{x},\hat{y},\hat{z}$$ Frame, AAE 532 Blake Lowe", 'Interpreter','Latex')
            xlabel('$$\hat x$$ direction [km]', 'Interpreter','Latex')
            ylabel('$$\hat y$$ direction [km]', 'Interpreter','Latex')
            zlabel('$$\hat z$$ direction [km]', 'Interpreter','Latex')
            axis equal

            % earth sphere
            if (draw_body)
                R_e = Body_Radius(obj.Body);
                [sph_x, sph_y, sph_z] = sphere(50);
                surf(sph_x*R_e, sph_y*R_e, sph_z*R_e)
                colormap summer
                shading interp
            end

            if (draw_pos)
                scatter3(obj.R_XYZ(1), obj.R_XYZ(2), obj.R_XYZ(3), 'black')
            end
            %xlim([-5e4, 5e4])
            %ylim([-5e4, 5e4])
            view([1,-1,1])
            hold off

        end

        function plot_XYZ_overlay(obj, n, orbit2, draw_body, draw_pos)
            ta_vec = linspace(0, 2*pi, n);
            r_mag_vec = obj.SLR./(1+obj.ECC*cos(ta_vec));
            r_e_vec = r_mag_vec.*cos(ta_vec);
            r_p_vec = r_mag_vec.*sin(ta_vec);
            r_eph_vec = [r_e_vec; r_p_vec; zeros(1,n)];
            r_xyz_vec = zeros(3, n);

            r2_mag_vec = orbit2.SLR./(1+orbit2.ECC*cos(ta_vec));
            r2_e_vec = r2_mag_vec.*cos(ta_vec);
            r2_p_vec = r2_mag_vec.*sin(ta_vec);
            r2_eph_vec = [r2_e_vec; r2_p_vec; zeros(1,n)];
            r2_xyz_vec = zeros(3, n);

            for i = 1:n
                r_xyz_vec(:,i) = eph2xyz(r_eph_vec(:,i), obj.AOP, obj.INC, obj.RAAN);
                r2_xyz_vec(:,i) = eph2xyz(r2_eph_vec(:,i), orbit2.AOP, orbit2.INC, orbit2.RAAN);
            end
            r_x_vec = r_xyz_vec(1,:);
            r_y_vec = r_xyz_vec(2,:);
            r_z_vec = r_xyz_vec(3,:);
            r2_x_vec = r2_xyz_vec(1,:);
            r2_y_vec = r2_xyz_vec(2,:);
            r2_z_vec = r2_xyz_vec(3,:);

            figure()
            hold on
            plot3(r_x_vec, r_y_vec, r_z_vec)
            plot3(r2_x_vec, r2_y_vec, r2_z_vec)
            title("Orbit Plot, $$\hat{x},\hat{y},\hat{z}$$ Frame, AAE 532 Blake Lowe", 'Interpreter','Latex')
            xlabel('$$\hat x$$ direction [km]', 'Interpreter','Latex')
            ylabel('$$\hat y$$ direction [km]', 'Interpreter','Latex')
            zlabel('$$\hat z$$ direction [km]', 'Interpreter','Latex')
            axis equal

            % earth sphere
            if (draw_body)
                R_e = Body_Radius(obj.Body);
                [sph_x, sph_y, sph_z] = sphere(50);
                surf(sph_x*R_e, sph_y*R_e, sph_z*R_e)
                colormap summer
                shading interp
            end

            if (draw_pos)
                scatter3(obj.R_XYZ(1), obj.R_XYZ(2), obj.R_XYZ(3), 'black')
            end
            %xlim([-5e4, 5e4])
            %ylim([-5e4, 5e4])
            view([1,-1,1])
            hold off

        end

        function plot_EPH_tikz(obj)
            % TODO
        end

        function plot_XYZ_tikz(obj)
            % TODO
        end
    end
end