function [out] = xyz2rth(in, theta, i, Omega)
% in, column vector in inertial frame to be rotated into orbit-fixed frame
% theta, argument of latitude (argument of periapsis + true anomaly
% i, inclination
% Omega, right ascension of the ascending node
R = [cos(Omega)*cos(theta) - sin(Omega)*cos(i)*sin(theta), sin(Omega)*cos(theta) + cos(Omega)*cos(i)*sin(theta), sin(i)*sin(theta);
    - cos(Omega)*sin(theta) - sin(Omega)*cos(i)*cos(theta), cos(Omega)*cos(i)*cos(theta) - sin(Omega)*sin(theta), cos(theta)*sin(i);
    sin(Omega)*sin(i), -cos(Omega)*sin(i), cos(i)];
out = R*in;
end
