function [mu] = Body_Grav_Param(body)
    % Returns gravitational parameter in  km^3/s^2
    switch body
        case 'Sun'
            mu = 132712440017.99;
        case 'Moon'
            mu = 4902.8005821478;
        case 'Mercury'
            mu = 22032.080486418;
        case 'Venus'
            mu = 324858.59882646;
        case 'Earth'
            mu = 398600.4415;
        case 'Mars'
            mu = 42828.314258067;
        case 'Jupiter'
            mu = 126712767.8578;
        case 'Saturn'
            mu = 37940626.061137;
        case 'Uranus'
            mu = 5794549.0070719;
        case 'Neptune'
            mu = 6836534.0638793;
        case 'Pluto'
            mu = 981.600887707;
        otherwise
            error("Body: '"+ body +"' does not have a defined Gravitational Parameter")
    end
end

