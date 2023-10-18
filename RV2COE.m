function [sma, ecc, argp, inc, raan, ta] = RV2COE(mu, r_, v_)
%RV2COE
% mu, gravitational parameter in km^2/s^2
% r_, position vector in inertial frame
% v_, velocity vector in inertial frame
r = norm(r_);
v = norm(v_);
h_ = cross(r_,v_);
h = norm(h_);
n_ = cross([0,0,1], h_);
n = norm(n_);

% eccentricity vector
e_ = ((v^2 - mu/r)*r_ - dot(r_, v_)*v_)/mu; 
e = norm(e_);
xi = v^2/2 - mu/r; % specific energy

% semimajor axis and semilatus rectum
if e ~= 1.0
    a = -mu/(2*xi);
    p = a*(1-e^2);
else
    a = inf;
    p = h^2/2;
end

% inclination
i = acos(h_(3)/h);

% right ascension of the ascending node
Omega = acos(n_(1)/n);
if n_(2) < 0
    Omega = 2*pi - Omega;
end

% argument of periapsis
omega = acos(dot(n_, e_)/(n*e));
if e_(3) < 0
    omega = 2*pi - omega;
end

% true anomaly
nu = acos(dot(e_, r_)/(e*r));
if dot(r_, v_) < 0
    nu = 2*pi - nu;
end


sma = a;
ecc = e;
inc = i;
argp = omega;
raan = Omega;
ta = nu;
end

