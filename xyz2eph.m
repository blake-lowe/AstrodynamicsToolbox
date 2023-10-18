function [out] = xyz2eph(in, omega, i, Omega)
% in, column vector in inertial frame to be rotated into orbit-fixed frame
% theta, argument of latitude (argument of periapsis + true anomaly
% i, inclination
% Omega, right ascension of the ascending node
R = [cos(Omega)*cos(omega) - sin(Omega)*cos(i)*sin(omega), sin(Omega)*cos(omega) + cos(Omega)*cos(i)*sin(omega), sin(i)*sin(omega);
    - cos(Omega)*sin(omega) - sin(Omega)*cos(i)*cos(omega), cos(Omega)*cos(i)*cos(omega) - sin(Omega)*sin(omega), cos(omega)*sin(i);
    sin(Omega)*sin(i), -cos(Omega)*sin(i), cos(i)];
out = R*in;
end
