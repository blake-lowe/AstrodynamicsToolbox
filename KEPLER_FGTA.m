function [outputArg1,outputArg2] = KEPLER_FGTA(inputArg1,inputArg2)
%KEPLER_FGTA Summary of this function goes here
%   Detailed explanation goes here
f_dot_alt = sqrt(mu/p)*tan((TA1-TA)/2)*((1-cos(TA1-TA))/p - 1/r1 - 1/r)
g_dot_alt = 1 - (r/p)*(1-cos(TA1 - TA));
TODO
end

