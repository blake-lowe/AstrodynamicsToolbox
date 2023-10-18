function [out] = rth2xyz(in, theta, i, Omega)
% in, column vector in orbit-fixed frame to be rotated into inertial frame
% argL, argument of latitude (argument of periapsis + true anomaly
% inc, inclination
% raan, right ascension of the ascending node
R = transpose([cos(Omega)*cos(theta) - sin(Omega)*cos(i)*sin(theta), sin(Omega)*cos(theta) + cos(Omega)*cos(i)*sin(theta), sin(i)*sin(theta);
    - cos(Omega)*sin(theta) - sin(Omega)*cos(i)*cos(theta), cos(Omega)*cos(i)*cos(theta) - sin(Omega)*sin(theta), cos(theta)*sin(i);
    sin(Omega)*sin(i), -cos(Omega)*sin(i), cos(i)]);
out = R*in;
end

