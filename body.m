classdef body
    %BODY Representation of a celestial body
    properties
        Name            % Name as a string
        Focus           % Focus as Body type
        RotPeriod       % Axial Rotational Period [rev/day]
        Radius          % Mean Equatorial Radius [km]
        Mu              % Gravitational Parameter [km^3/s^2]
        SMA             % Semi-major Axis of Orbit [km]
        Period          % Orbital Period [s]
        ECC             % Eccentricity of Orbit [1]
        INC             % Inclination of Orbit [rad]
    end
    
    methods
        function obj = body(varargin)
            % body(name)
            % body(name, focus, rotPeriod, radius, mu, sma, period, ecc, inc)
            if length(varargin) == 1 % From name
                obj = body.bodyFromName(varargin{1});
            elseif length(varargin) == 9 % Fully defined
                obj.Name = varargin{1};
                obj.Focus = varargin{2};
                obj.RotPeriod = varargin{3};
                obj.Radius = varargin{4};
                obj.Mu = varargin{5};
                obj.SMA = varargin{6};
                obj.Period = varargin{7};
                obj.ECC = varargin{8};
                obj.INC = varargin{9};
            else
                error('Incorrect number of inputs')
            end
        end
    end

    methods (Static)
        
        function obj = bodyFromName(name)
            % Construct a body object from reference data.
            % Data from JPLJs ephemerides file de405.spk and
            % https://ssd.jpl.nasa.gov/planets/approx_pos.html
            switch name
                case 'Sun'
                    focus = NaN;
                    rotPeriod = 0.0394011;
                    radius = 695990;
                    mu = 132712440017.99;
                    sma = NaN;
                    period = NaN;
                    ecc = NaN;
                    inc = NaN;
                case 'Moon'
                    focus = body('Earth');
                    rotPeriod = 0.0366004;
                    radius = 1738.2;
                    mu = 4902.8005821478;
                    sma = 384400; % around Earth
                    period = 2360592;
                    ecc = 0.0554;
                    inc = deg2rad(5.16);
                case 'Mercury'
                    focus = body('Sun');
                    rotPeriod = 0.0170514;
                    radius = 2439.7;
                    mu = 22032.080486418;
                    sma = 57909101;
                    period = 7600537;
                    ecc = 0.20563661;
                    inc = deg2rad(7.00497902);
                case 'Venus'
                    focus = body('Sun');
                    rotPeriod = -0.0041149;
                    radius = 6051.9;
                    mu = 324858.59882646;
                    sma = 108207284;
                    period = 19413722;
                    ecc = 0.00676399;
                    inc = deg2rad(3.39465605);
                case 'Earth'
                    focus = body('Sun');
                    rotPeriod = 1.0027378;
                    radius = 6378.1363;
                    mu = 398600.4415;
                    sma = 149597898;
                    period = 31558205;
                    ecc = 0.01673163;
                    inc = deg2rad(0.00001531);
                case 'Mars'
                    focus = body('Sun');
                    rotPeriod = 0.9747000;
                    radius = 3397;
                    mu = 42828.314258067;
                    sma = 227944135;
                    period = 59356281;
                    ecc = 0.09336511;
                    inc = deg2rad(1.84969142);
                case 'Jupiter'
                    focus = body('Sun');
                    rotPeriod = 2.4181573;
                    radius = 71492;
                    mu = 126712767.8578;
                    sma = 778279959;
                    period = 374479305;
                    ecc = 0.04853590;
                    inc = deg2rad(1.30439695);
                case 'Saturn'
                    focus = body('Sun');
                    rotPeriod = 2.2522053;
                    radius = 60268;
                    mu = 37940626.061137;
                    sma = 1427387908;
                    period = 930115906;
                    ecc = 0.05550825;
                    inc = deg2rad(2.48599187);
                case 'Uranus'
                    focus = body('Sun');
                    rotPeriod = -1.3921114;
                    radius = 25559;
                    mu = 5794549.0070719;
                    sma = 2870480873;
                    period = 2652503938;
                    ecc = 0.04685740;
                    inc = deg2rad(0.77263783);
                case 'Neptune'
                    focus = body('Sun');
                    rotPeriod = 1.4897579;
                    radius = 25269;
                    mu = 6836534.0638793;
                    sma = 4498337290;
                    period = 5203578080;
                    ecc = 0.00895439;
                    inc = deg2rad(1.77004347);
                case 'Pluto'
                    focus = body('Sun');
                    rotPeriod = -0.1565620;
                    radius = 1162;
                    mu = 981.600887707;
                    sma = 5907150229;
                    period = 7830528509;
                    ecc = 0.24885238;
                    inc = deg2rad(17.14001206);
                otherwise
                    error("Body: '"+ body +"' properties have not been defined.")
            end
            obj = body(name, focus, rotPeriod, radius, mu, sma, period, ecc, inc);
        end
    end
end
            
