function [r_orb, v_orb, r_rot, v_rot] = RVTH2COMPONENTS(r, v, FPA, TA)
    r_rot = [r, 0, 0]';
    v_rot = [v*sin(FPA), v*cos(FPA), 0]';
    
    % Rotation matrix from orbit frame to rotating frame
    R = [cos(TA) sin(TA) 0;
        -sin(TA) cos(TA) 0;
        0 0 1];

    r_orb = R'*r_rot;
    v_orb = R'*v_rot;
end

