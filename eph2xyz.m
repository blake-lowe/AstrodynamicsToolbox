function [out] = eph2xyz(in, omega, i, Omega)
% in, column vector in orbit-fixed frame to be rotated into inertial frame
% argL, argument of latitude (argument of periapsis + true anomaly
% inc, inclination
% raan, right ascension of the ascending node
R = transpose([cos(Omega)*cos(omega) - sin(Omega)*cos(i)*sin(omega), sin(Omega)*cos(omega) + cos(Omega)*cos(i)*sin(omega), sin(i)*sin(omega);
    - cos(Omega)*sin(omega) - sin(Omega)*cos(i)*cos(omega), cos(Omega)*cos(i)*cos(omega) - sin(Omega)*sin(omega), cos(omega)*sin(i);
    sin(Omega)*sin(i), -cos(Omega)*sin(i), cos(i)]);
out = R*in;
end

