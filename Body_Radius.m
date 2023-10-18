function [radius] = Body_Radius(body)
    % Returns Mean Equatorial Radius in km
    switch body
        case 'Sun'
            radius = 695990;
        case 'Moon'
            radius = 1738.2;
        case 'Mercury'
            radius = 2439.7;
        case 'Venus'
            radius = 6051.9;
        case 'Earth'
            radius = 6378.1363;
        case 'Mars'
            radius = 3397;
        case 'Jupiter'
            radius = 71492;
        case 'Saturn'
            radius = 60268;
        case 'Uranus'
            radius = 25559;
        case 'Neptune'
            radius = 25269;
        case 'Pluto'
            radius = 1162;
        otherwise
            error("Body: '"+ body +"' does not have a defined Radius")
    end
end

