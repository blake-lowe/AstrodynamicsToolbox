function [out] = rth2eph(in, TA)
    % Rotation matrix from orbit frame to rotating frame
    R = [cos(TA) sin(TA) 0;
        -sin(TA) cos(TA) 0;
        0 0 1];

    out = R'*in;
end

