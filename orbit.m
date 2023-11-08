classdef orbit
    %ORBIT Representation of a spacecraft on a conic orbit
    
    properties % Distances in km, Angles in rad, Time in s
        Body        % Primary body, as a body class
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
        R_DOT       % Range Rate
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
            %orbit('COE', Body, SMA, ECC, AOP, INC, RAAN, TA)
            %orbit('RARP', Body, RA, RP, AOP, INC, RAAN, TA)
            %orbit('RV', Body, R_XYZ, V_XYZ)
            %orbit('RRP', Body, R1_XYZ, R2_XYZ, P) not implemented
            %orbit('RRT', Body, R1_XYZ, R2_XYZ, T21) not implemented
            if ischar(varargin{2})
                obj.Body = body(varargin{2});
            else
                obj.Body = varargin{2};
            end

            if strcmp(varargin{1},'COE') % Keplerian Orbital Elements
                obj.SMA = varargin{3};
                obj.ECC = varargin{4};
                obj.AOP = varargin{5};
                obj.INC = varargin{6};
                obj.RAAN = varargin{7};
                obj.TA = varargin{8};
            elseif strcmp(varargin{1},'RARP') % Apoapsis and Periapsis
                RA = varargin{3};
                RP = varargin{4};
                obj.SMA = (RA + RP)/2;
                obj.ECC = (RA - RP)/(RA + RP);
                obj.AOP = varargin{5};
                obj.INC = varargin{6};
                obj.RAAN = varargin{7};
                obj.TA = varargin{8};
            elseif strcmp(varargin{1}, 'RV') % Position and velocity
                [obj.SMA, obj.ECC, obj.AOP, obj.INC, obj.RAAN, obj.TA] =...
                    RV2COE(obj.Body.Mu, varargin{3}, varargin{4});
            elseif strcmp(varargin{1}, 'RRP') % Two Positions and P
                % TODO using f and g functions to get v1 then RV2COE(r1,v1)

            elseif strcmp(varargin{1}, 'Body') % todo date and ephemerides
                body = varargin{2};
                obj = orbit('COE', body.Focus, body.SMA, body.ECC, 0, body.INC, 0, 0);
            else
                error('Invalid input mode specified')
            end

            obj = obj.update_properties();
        end

        function obj = propagate_TA(obj, dTA)
            if (dTA < 0 && obj.Type ~= 0)
                warning('Backwards propagation occurred on a non-elliptical orbit')
            end
            obj.TA = obj.TA + dTA;
            TP_before = obj.TP;
            obj = obj.update_properties();
            TP_after = obj.TP;
            if obj.Type == 0
                propagationTime = mod(TP_after-TP_before, obj.Period)
            else
                propagationTime = TP_after-TP_before
            end
        end

        function obj = propagate_toTA(obj, TA)
            if (obj.Type == 0) && (TA > 2*pi || TA < 0)
                error('Invalid New TA for elliptical orbit. Should be [0,2*pi]')
            end
            if (obj.Type == 1) && (TA > pi || TA < -pi)
                error('Invalid New TA for parabolic orbit. Should be [-pi,pi]')
            end
            %TODO error case for hyperbolic orbit (use ta inf)

            % Wrap TA for elliptical orbits
            if obj.ECC < 1 && obj.TA > TA
                TA = TA + 2*pi;
            end

            dTA = TA - obj.TA;
            obj = obj.propagate_TA(dTA);
        end

        function obj = propagate_toR(obj, R, Ascending)
            %propagate_toR(R, Ascending)
            % Propagate the orbit a given Radius
            % R, radius
            % Ascending, logical true if ascending
            if ~isnan(obj.RA) && R > obj.RA
                error('Radius is too high, it is never reached')
            end
            if (R < obj.RP)
                error('Radius is too low, it is never reached')
            end

            TA = acos((1/obj.ECC)*(obj.SLR/R - 1));
            if (~Ascending)
                TA = 2*pi - TA;
            end
            dTA = TA - obj.TA;
            obj = obj.propagate_TA(dTA);

        end

        function obj = propagate_EA(obj, dEA)
            % Calculate dt
            M0 = obj.MA;
            M1 = Anomaly.EA2MA(obj.EA + dEA, obj.ECC);
            dt = (M1 - M0)/obj.MM;

            % Use f and g functions to get new R,V
            [R1_XYZ, V1_XYZ] = KEPLER_FGEA(obj.Body.Mu, obj.SMA, dt, dEA, obj.R_XYZ, obj.V_XYZ);

            % Construct the new orbit RV constructor
            TP_before = obj.TP;
            obj = orbit('RV', obj.Body, R1_XYZ, V1_XYZ);
            TP_after = obj.TP;
            propagationTime = mod(TP_after-TP_before, obj.Period)
        end

        function obj = propagate_toEA(obj, EA)
            if (obj.Type == 0) && (EA > 2*pi || EA < 0)
                error('Invalid New EA for elliptical orbit. Should be [0,2*pi]')
            end
            if (obj.Type ~= 0)
                error('Cannot propagate by EA for non-elliptical orbit')
            end

            % Wrap TA for elliptical orbits
            if obj.EA > EA
                EA = EA + 2*pi;
            end

            dEA = EA - obj.EA;
            obj = obj.propagate_EA(dEA);
        end


        function obj = propagate_Time(obj, dt)
            obj.XA = obj.XA + dt*obj.MM;
            if obj.Type == 0 % Ellipse
                obj.MA = obj.XA;
                obj.EA = Anomaly.MA2EA(obj.MA, obj.ECC, 1e-12);
                obj.TA = mod(2*atan(sqrt((1+obj.ECC)/(1-obj.ECC))*tan(obj.EA/2)), 2*pi);
            elseif obj.Type == 1 % Parabola
                obj.OA = obj.XA;
                obj.BA = Anomaly.OA2BA(obj.OA);
                obj.TA = 2*atan(obj.BA);
            else % Hyperbola
                obj.NA = obj.XA;
                obj.HA = Anomaly.NA2HA(obj.NA, obj.ECC, 1e-12);
                obj.TA = mod(2*atan(sqrt((1+obj.ECC)/(1-obj.ECC))*tanh(obj.HA/2)), 2*pi);
            end

            obj = obj.update_properties();
        end

        function obj = maneuver_execute(obj, maneuver1)
            obj = obj.maneuver_XYZ(maneuver1.dV_XYZ);
        end

        function obj = maneuver_XYZ(obj, dV_XYZ)
            if length(dV_XYZ) ~= 3
                error('dV must have length = 3')
            end
            dV = norm(dV_XYZ)
            obj = orbit('RV', obj.Body, obj.R_XYZ, obj.V_XYZ + dV_XYZ);
        end

        function obj = update_properties(obj)
            %UPDATE_PROPERTIES, calculates all derived non COE properties

            %MU          % Gravitational parameter of Body
            obj.MU = obj.Body.Mu;

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
                obj.HA = mod(2*atanh(sqrt((obj.ECC-1)/(obj.ECC+1))*tan(obj.TA/2)), 2*pi);
                % NA          % Hyperbolic Mean Anomaly
                obj.NA = obj.ECC*sinh(obj.HA) - obj.HA;
                obj.MM = sqrt(obj.MU/(-obj.SMA^3));
                obj.XA = obj.NA;
                obj.RP = obj.SMA*(1-obj.ECC);
                obj.RA = NaN;
            end

            if obj.RP < obj.Body.Radius
                warning('Orbit intersects body surface.')
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
            obj.FPA = norm(acos(obj.H/(obj.R*obj.V)));
            if obj.Ascending
                obj.FPA = -obj.FPA;
            end

            % R_RTH       % Position vector in rotating orbit frame
            obj.R_RTH = [obj.R; 0; 0];
            % V_RTH       % Velocity vector in rotating orbit frame
            obj.V_RTH = [obj.V*sin(obj.FPA); obj.V*cos(obj.FPA); 0];

            % R_DOT       % Range Rate, projection of V_RTH onto R_RTH
            obj.R_DOT = norm(dot(obj.R_RTH, obj.V_RTH)/obj.R);

            % R_EPH       % Position vector in fixed orbit frame
            obj.R_EPH = Frame.rth2eph(obj.R_RTH, obj.TA);
            % V_EPH       % Velcotiy vector in fixed orbit frame
            obj.V_EPH = Frame.rth2eph(obj.V_RTH, obj.TA);
            
            % R_XYZ       % Position vector in inertial frame
            obj.R_XYZ = Frame.rth2xyz(obj.R_RTH, obj.AOL, obj.INC, obj.RAAN);
            % V_XYZ       % Velocity vector in inertial frame
            obj.V_XYZ = Frame.rth2xyz(obj.V_RTH, obj.AOL, obj.INC, obj.RAAN);
            if isnan(obj.RAAN) % in an equatorial orbit
                obj.R_XYZ = Frame.rth2xyz(obj.R_RTH, obj.AOL, obj.INC, 0);
                obj.V_XYZ = Frame.rth2xyz(obj.V_RTH, obj.AOL, obj.INC, 0);
            end

            % H_XYZ       % Angular Momentum Vector in inertial frame
            obj.H_XYZ = cross(obj.R_XYZ, obj.V_XYZ);

            % E_XYZ       % Eccentricity Vector in inertial frame
            obj.E_XYZ = (1/obj.MU)*((obj.V^2 - obj.MU/obj.R)*obj.R_XYZ - dot(obj.R_XYZ, obj.V_XYZ)*obj.V_XYZ);

            % N_XYZ       % Nodal Vector in inertial frame
            obj.N_XYZ = cross([0;0;1], obj.H_XYZ);
            obj.N_XYZ = obj.N_XYZ./norm(obj.N_XYZ);
        end

        function plot_EPH(obj, n, draw_body, draw_pos, draw_vecs, draw_apsides, draw_nodes, draw_aux, legend_labels) % Only works for ellipses for now
            if obj.ECC >= 1
                ta_lim = min(acos(-1/obj.ECC), max(2*pi/3, abs(obj.TA*1.01)));
                ta_vec = linspace(-ta_lim, ta_lim, n);
            else
                ta_vec = linspace(0, 2*pi, n);
            end
            
            r_mag_vec = obj.SLR./(1+obj.ECC*cos(ta_vec));
            r_e_vec = r_mag_vec.*cos(ta_vec);
            r_p_vec = r_mag_vec.*sin(ta_vec);
            aux_e_vec = -obj.SMA*obj.ECC + obj.SMA*cos(ta_vec);
            aux_p_vec = obj.SMA*sin(ta_vec);
            
            figure()
            hold on
            plot(r_e_vec, r_p_vec,'LineWidth',1.5)
            xlabel("$$\hat{e}$$ direction [km]",'Interpreter','Latex')
            ylabel("$$\hat{p}$$ direction [km]",'Interpreter','Latex')
            title("Orbit Plot, $$\hat{e},\hat{p},\hat{h}$$ Frame, AAE 532 Blake Lowe", 'Interpreter','Latex')
            axis equal
            if (draw_aux)
                plot(aux_e_vec, aux_p_vec)
            end
            if (draw_body)
                p = nsidedpoly(1000, 'Center', [0 0], 'Radius', obj.Body.Radius);
                plot(p, 'FaceColor', 'b')
            end
            scatter(0, 0, 'black', 'x') % focus
            scatter(-obj.SMA*obj.ECC, 0, 'black', '+') %center
            if (draw_apsides)
                plot([-obj.RA,obj.RP], [0,0], 'black--')
            end
            if (draw_nodes)
                % TODO
                %N_EPH = Frame.xyz2eph(N_XYZ, obj.AOP, obj.INC, obj.RAAN);
                plot([-cos(-obj.AOP), cos(-obj.AOP)]*obj.SMA, [-sin(-obj.AOP), sin(-obj.AOP)]*3*obj.SMA, 'red--')
            end
            if (draw_pos)
                scatter(obj.R_EPH(1), obj.R_EPH(2), 'black')
            end
            if (draw_vecs)
                % TODO unit vectors, vectors, and angles
                plot([0,obj.R_EPH(1)], [0,obj.R_EPH(2)], 'black')
            end
            %ylim([min(r_p_vec) - 0.1*obj.SMA, max(r_p_vec) + 0.1*obj.SMA])

            if exist('legend_labels', 'var')
                legend(legend_labels)
            end

            hold off
        end

        function plot_EPH_overlay(obj, n, orbit2, draw_body, draw_pos, draw_vecs, draw_apsides, legend_labels)
            if obj.ECC >= 1
                ta_lim = min(acos(-1/obj.ECC), max(2*pi/3, abs(obj.TA*1.01)));
                ta_vec = linspace(-ta_lim, ta_lim, n);
            else
                ta_vec = linspace(0, 2*pi, n);
            end
            r_mag_vec = obj.SLR./(1+obj.ECC*cos(ta_vec));
            r_e_vec = r_mag_vec.*cos(ta_vec);
            r_p_vec = r_mag_vec.*sin(ta_vec);

            if orbit2.ECC >= 1
                ta_lim = min(acos(-1/orbit2.ECC), max(2*pi/3, abs(orbit2.TA*1.01)));
                ta_vec = linspace(-ta_lim, ta_lim, n);
            else
                ta_vec = linspace(0, 2*pi, n);
            end
            r2_mag_vec = orbit2.SLR./(1+orbit2.ECC*cos(ta_vec));
            r2_e_vec = zeros(1, n);
            r2_p_vec = zeros(1, n);

            for i = 1:n
                r2_xyz = Frame.rth2xyz([r2_mag_vec(i);0;0], orbit2.AOP + ta_vec(i), orbit2.INC, orbit2.RAAN);
                r2_eph1 = Frame.xyz2eph(r2_xyz, obj.AOP, obj.INC, obj.RAAN);
                r2_e_vec(i) = r2_eph1(1);
                r2_p_vec(i) = r2_eph1(2);
            end


            figure()
            hold on
            plot(r_e_vec, r_p_vec,'LineWidth',1.5)
            plot(r2_e_vec, r2_p_vec,'LineWidth',1.5)
            xlabel("$$\hat{e}$$ direction [km]",'Interpreter','Latex')
            ylabel("$$\hat{p}$$ direction [km]",'Interpreter','Latex')
            title("Overlaid Orbit Plot, $$\hat{e},\hat{p},\hat{h}$$ Frame, AAE 532 Blake Lowe", 'Interpreter','Latex')
            axis equal
            if (draw_body)
                p = nsidedpoly(1000, 'Center', [0 0], 'Radius', obj.Body.Radius);
                plot(p, 'FaceColor', 'b')
            end
            scatter(0, 0, 'black', 'x') % focus
            scatter(-obj.SMA*obj.ECC, 0, 'black', '+') %center
            if (draw_apsides)
                plot([-obj.RA,obj.RP], [0,0], 'black--')
                peri_xyz  = Frame.eph2xyz([orbit2.RP;0;0], orbit2.AOP + 0, orbit2.INC, orbit2.RAAN);
                peri_eph1 = Frame.xyz2eph(peri_xyz, obj.AOP, obj.INC, obj.RAAN);
                apo_xyz  = Frame.eph2xyz([-orbit2.RA;0;0], orbit2.AOP + 0, orbit2.INC, orbit2.RAAN);
                apo_eph1 = Frame.xyz2eph(apo_xyz, obj.AOP, obj.INC, obj.RAAN);
                plot([peri_eph1(1),apo_eph1(1)], [peri_eph1(2),apo_eph1(2)], 'black--')
            end
            if (draw_pos)
                scatter(obj.R_EPH(1), obj.R_EPH(2), 'black')
                R2_XYZ = Frame.rth2xyz([orbit2.R;0;0], orbit2.AOL, orbit2.INC, orbit2.RAAN);
                R2_EPH1 = Frame.xyz2eph(R2_XYZ, obj.AOP, obj.INC, obj.RAAN);
                scatter(R2_EPH1(1), R2_EPH1(2), 'black')
            end
            if (draw_vecs)
                % TODO unit vectors, vectors, and angles
                plot([0,obj.R_EPH(1)], [0,obj.R_EPH(2)], 'black')
                R2_XYZ = Frame.rth2xyz([orbit2.R;0;0], orbit2.AOL, orbit2.INC, orbit2.RAAN);
                R2_EPH1 = Frame.xyz2eph(R2_XYZ, obj.AOP, obj.INC, obj.RAAN);
                plot([0,R2_EPH1(1)], [0,R2_EPH1(2)], 'black')
            end

            minx = min([r_e_vec, r2_e_vec]);
            maxx = max([r_e_vec, r2_e_vec]);
            miny = min([r_p_vec, r2_p_vec]);
            maxy = max([r_p_vec, r2_p_vec]);
            biga = max(obj.SMA, orbit2.SMA);

            xlim([minx, maxx]+[-0.1, 0.1]*biga)
            ylim([miny, maxy]+[-0.1, 0.1]*biga)
            axis equal

            if exist('legend_labels', 'var')
                legend(legend_labels)
            end

            hold off
        end

        function plot_XYZ(obj, n, draw_body, draw_plane, draw_pos, draw_apsides, legend_labels)
            if obj.ECC >= 1
                ta_lim = min(acos(-1/obj.ECC), max(2*pi/3, abs(obj.TA*1.01)));
                ta_vec = linspace(-ta_lim, ta_lim, n);
            else
                ta_vec = linspace(0, 2*pi, n);
            end

            r_mag_vec = obj.SLR./(1+obj.ECC*cos(ta_vec));
            r_e_vec = r_mag_vec.*cos(ta_vec);
            r_p_vec = r_mag_vec.*sin(ta_vec);
            r_eph_vec = [r_e_vec; r_p_vec; zeros(1,n)];
            r_xyz_vec = zeros(3, n);
            for i = 1:n
                r_xyz_vec(:,i) = Frame.eph2xyz(r_eph_vec(:,i), obj.AOP, obj.INC, obj.RAAN);
            end
            r_x_vec = r_xyz_vec(1,:);
            r_y_vec = r_xyz_vec(2,:);
            r_z_vec = r_xyz_vec(3,:);

            figure()
            hold on
            plot3(r_x_vec, r_y_vec, r_z_vec,'LineWidth',1.5)
            title("Orbit Plot, $$\hat{x},\hat{y},\hat{z}$$ Frame, AAE 532 Blake Lowe", 'Interpreter','Latex')
            xlabel('$$\hat x$$ direction [km]', 'Interpreter','Latex')
            ylabel('$$\hat y$$ direction [km]', 'Interpreter','Latex')
            zlabel('$$\hat z$$ direction [km]', 'Interpreter','Latex')
            axis equal

            % earth sphere
            if (draw_body)
                R_e = obj.Body.Radius;
                [sph_x, sph_y, sph_z] = sphere(50);
                surf(sph_x*R_e, sph_y*R_e, sph_z*R_e)
                colormap summer
                shading interp
            end

            if (draw_pos)
                scatter3(obj.R_XYZ(1), obj.R_XYZ(2), obj.R_XYZ(3), 'black', 'filled')
            end

            if (draw_apsides)
                peri1_xyz = Frame.eph2xyz([obj.RP;0;0], obj.AOP, obj.INC, obj.RAAN);
                apo1_xyz = Frame.eph2xyz([-obj.RA;0;0], obj.AOP, obj.INC, obj.RAAN);
                plot3([peri1_xyz(1),apo1_xyz(1)], [peri1_xyz(2),apo1_xyz(2)], [peri1_xyz(3), apo1_xyz(3)], 'black--')
            end
            
            buf = 0.1*obj.SMA;
            minx = min([r_x_vec, -obj.Body.Radius])-buf;
            maxx = max([r_x_vec, obj.Body.Radius])+buf;
            miny = min([r_y_vec, -obj.Body.Radius])-buf;
            maxy = max([r_y_vec, obj.Body.Radius])+buf;
            minz = min([r_z_vec, -obj.Body.Radius])-buf;
            maxz = max([r_z_vec, obj.Body.Radius])+buf;

            % fundamental plane
            if (draw_plane)
                p = patch([minx, maxx, maxx, minx], [miny, miny, maxy, maxy], [0,0,0,0], 'black');
                set(p,'facealpha',0.5);
            end

            xlim([minx, maxx])
            ylim([miny, maxy])
            zlim([minz, maxz])

            view([1,-1,1])
            hold off

            if exist('legend_labels', 'var')
                legend(legend_labels)
            end
        end

        function plot_XYZ_overlay(obj, n, orbits, draw_body, draw_plane, draw_pos, draw_apsides, legend_labels)
            if obj.ECC >= 1
                ta_lim = min(acos(-1/obj.ECC), max(2*pi/3, abs(obj.TA*1.01)));
                ta_vec = linspace(-ta_lim, ta_lim, n);
            else
                ta_vec = linspace(0, 2*pi, n);
            end
            r_mag_vec = obj.SLR./(1+obj.ECC*cos(ta_vec));
            r_e_vec = r_mag_vec.*cos(ta_vec);
            r_p_vec = r_mag_vec.*sin(ta_vec);
            r_eph_vec = [r_e_vec; r_p_vec; zeros(1,n)];
            r_xyz_vec = zeros(3, n);

            for i = 1:n
                r_xyz_vec(:,i) = Frame.eph2xyz(r_eph_vec(:,i), obj.AOP, obj.INC, obj.RAAN);
            end

            r_x_vec = r_xyz_vec(1,:);
            r_y_vec = r_xyz_vec(2,:);
            r_z_vec = r_xyz_vec(3,:);

            minx = min([r_x_vec, -obj.Body.Radius]);
            maxx = max([r_x_vec, obj.Body.Radius]);
            miny = min([r_y_vec, -obj.Body.Radius]);
            maxy = max([r_y_vec, obj.Body.Radius]);
            minz = min([r_z_vec, -obj.Body.Radius]);
            maxz = max([r_z_vec, obj.Body.Radius]);

            figure()
            hold on
            plot3(r_x_vec, r_y_vec, r_z_vec,'LineWidth',1.5)
            title("Orbit Plot, $$\hat{x},\hat{y},\hat{z}$$ Frame, AAE 532 Blake Lowe", 'Interpreter','Latex')
            xlabel('$$\hat x$$ direction [km]', 'Interpreter','Latex')
            ylabel('$$\hat y$$ direction [km]', 'Interpreter','Latex')
            zlabel('$$\hat z$$ direction [km]', 'Interpreter','Latex')
            axis equal

            biga = obj.SMA;

            for o = 1:length(orbits)
                if orbits(o).ECC >= 1
                    ta_lim = min(acos(-1/orbits(o).ECC), max(2*pi/3, abs(orbits(o).TA*1.01)));
                    ta_vec = linspace(-ta_lim, ta_lim, n);
                else
                    ta_vec = linspace(0, 2*pi, n);
                end
                ro_mag_vec = orbits(o).SLR./(1+orbits(o).ECC*cos(ta_vec));
                ro_e_vec = ro_mag_vec.*cos(ta_vec);
                ro_p_vec = ro_mag_vec.*sin(ta_vec);
                ro_eph_vec = [ro_e_vec; ro_p_vec; zeros(1,n)];
                ro_xyz_vec = zeros(3, n);

                for i = 1:n
                    ro_xyz_vec(:,i) = Frame.eph2xyz(ro_eph_vec(:,i), orbits(o).AOP, orbits(o).INC, orbits(o).RAAN);
                end

                ro_x_vec = ro_xyz_vec(1,:);
                ro_y_vec = ro_xyz_vec(2,:);
                ro_z_vec = ro_xyz_vec(3,:);

                plot3(ro_x_vec, ro_y_vec, ro_z_vec,'LineWidth',1.5)
                
                if orbits(o).SMA > biga
                    biga = orbits(o).SMA;
                end
                if min(ro_x_vec) < minx
                    minx = min(ro_x_vec);
                end
                if max(ro_x_vec) > maxx
                    maxx = max(ro_x_vec);
                end
                if min(ro_y_vec) < miny
                    miny = min(ro_y_vec);
                end
                if max(ro_y_vec) > maxy
                    maxy = max(ro_y_vec);
                end
                if min(ro_z_vec) < minz
                    minz = min(ro_z_vec);
                end
                if max(ro_z_vec) > maxz
                    maxz = max(ro_z_vec);
                end

            end

            % body sphere
            if (draw_body)
                R_e = obj.Body.Radius;
                [sph_x, sph_y, sph_z] = sphere(50);
                v = surf(sph_x*R_e, sph_y*R_e, sph_z*R_e);
                colormap summer
                shading interp
                set(v,'facealpha',0.65);
            end

            if (draw_pos)
                scatter3(obj.R_XYZ(1), obj.R_XYZ(2), obj.R_XYZ(3), 'black', 'filled')

                for o = 1:length(orbits)
                    Ro_XYZ = Frame.rth2xyz([orbits(o).R;0;0], orbits(o).AOL, orbits(o).INC, orbits(o).RAAN);
                    scatter3(Ro_XYZ(1), Ro_XYZ(2), Ro_XYZ(3), 'black','filled')
                end
            end

            if (draw_apsides)
                peri1_xyz = Frame.eph2xyz([obj.RP;0;0], obj.AOP, obj.INC, obj.RAAN);
                apo1_xyz = Frame.eph2xyz([-obj.RA;0;0], obj.AOP, obj.INC, obj.RAAN);
                plot3([peri1_xyz(1),apo1_xyz(1)], [peri1_xyz(2),apo1_xyz(2)], [peri1_xyz(3), apo1_xyz(3)], 'black--')

                for o = 1:length(orbits)
                    perio_xyz = Frame.eph2xyz([orbits(o).RP;0;0], orbits(o).AOP, orbits(o).INC, orbits(o).RAAN);
                    apoo_xyz  = Frame.eph2xyz([-orbits(o).RA;0;0], orbits(o).AOP, orbits(o).INC, orbits(o).RAAN);
                    plot3([perio_xyz(1),apoo_xyz(1)], [perio_xyz(2),apoo_xyz(2)], [perio_xyz(3), apoo_xyz(3)], 'black--')
                end
            end
            
            buf = 0.1*biga;
            
            xbounds = [minx-buf, maxx+buf];
            ybounds = [miny-buf, maxy+buf];
            zbounds = [minz-buf, maxz+buf];

            % fundamental plane
            if (draw_plane)
                p = patch([xbounds(1), xbounds(2), xbounds(2), xbounds(1)], [ybounds(1), ybounds(1), ybounds(2), ybounds(2)], [0,0,0,0], 'black');
                set(p,'facealpha',0.5);
            end

            xlim(xbounds)
            ylim(ybounds)
            zlim(zbounds)

            view([1,-1,1])
            hold off

            if exist('legend_labels', 'var')
                legend(legend_labels)
            end
        end

        function plot_EPH_tikz(obj, draw_body, draw_pos, draw_apsides, draw_vectors, draw_angles)
            % plot_EPH_tikz(obj, draw_body, draw_pos, draw_apsides, draw_vectors, draw_angles)

            v_s = 3;% higher is smaller velocity vectors
            r_s = 1e4;
            if(obj.Type ~= 0)
                error('Non-elliptical orbits not supported yet')
            end
            fprintf('\\begin{center}\n')
            fprintf('\\begin{tikzpicture}[scale = 1]\n')
            fprintf('\\coordinate (focus) at (0,0);\n')
            fprintf('\\coordinate (sc) at (%0.4f,%0.4f);\n', obj.R_EPH(1)/r_s, obj.R_EPH(2)/r_s)
            b = obj.SMA*sqrt(1-obj.ECC^2);
            fprintf('\\draw [line width = 1.2pt, blue] (%0.4f,0) ellipse (%0.4f and %0.4f);\n', -obj.SMA*obj.ECC/r_s, obj.SMA/r_s, b/r_s)
            if (draw_body)
                fprintf('\\draw [fill=white] (focus) circle [radius = %0.4f];\n', obj.Body.Radius/r_s)
            end
            if (draw_apsides)
                fprintf('\\draw [dashed] (%0.4f,0) -- (%0.4f,0);\n', -obj.RA/1e4, obj.RP/r_s)
            end
            fprintf('\\node at (focus) [shift={(0,0)}] {$\\%s$};\n',obj.Body.Name)
            e_hat = [obj.SMA/3e4;0;0];
            p_hat = [0;obj.SMA/(3*r_s);0];
            r_hat = Frame.rth2eph(e_hat, obj.TA);
            t_hat = Frame.rth2eph(p_hat, obj.TA);
            if (draw_vectors)
                fprintf('\\draw [teal, thick, ->](focus) -- node [midway, above right] {$\\bar{r}$} (sc);\n')
                fprintf('\\draw [teal, ->](sc) -- ($(sc) + (%0.4f,%0.4f)$) node [above] {$\\bar{v}$};\n',obj.V_EPH(1)/v_s, obj.V_EPH(2)/v_s)
                fprintf('\\draw [red, ->](focus) -- ($(focus) + (%0.4f,%0.4f)$) node [below right] {$\\hat{e}$};\n', e_hat(1), e_hat(2))
                fprintf('\\draw [red, ->](focus) -- ($(focus) + (%0.4f,%0.4f)$) node [left] {$\\hat{p}$};\n', p_hat(1), p_hat(2))
                fprintf('\\draw [red, ->](sc) -- ($(sc) + (%0.4f,%0.4f)$) node [above right] {$\\hat{r}$};\n', r_hat(1), r_hat(2))
                fprintf('\\draw [red, ->](sc) -- ($(sc) + (%0.4f,%0.4f)$) node [below] {$\\hat{\\theta}$};\n', t_hat(1), t_hat(2))
                fprintf('\\draw [gray, dashed](sc) -- ($(sc) - (%0.4f,%0.4f)$) node [above] {l.h.};\n', t_hat(1), t_hat(2))
            end
            if (draw_angles)
                fprintf('\\draw[ultra thick,orange, |->] ([shift=(0:%0.4f)]focus) arc (0:%0.4f:%0.4f) node [midway, above right]{$\\theta^*$};\n', obj.RP/(1.5*r_s), rad2deg(obj.TA), obj.RP/(1.5*r_s))
                ang1 = rad2deg(obj.TA+pi/2);
                ang2 = ang1 - rad2deg(obj.FPA);
                fprintf('\\draw[ultra thick,orange, |->] ([shift=(%0.4f:%0.4f)]sc) arc (%0.4f:%0.4f:%0.4f) node [midway, below left]{$\\gamma$};\n', ang1, obj.V/(2.5*v_s), ang1, ang2, obj.V/(2.5*v_s))
            end
            if (draw_pos)
                fprintf('\\draw (sc) circle [radius = 0.1];%% node [right] {$s/c$};\n')
            end
            fprintf('\\end{tikzpicture}\n')
            fprintf('\\end{center}\n')
        end

        function plot_XYZ_tikz(obj)
            % TODO
        end
    end
end