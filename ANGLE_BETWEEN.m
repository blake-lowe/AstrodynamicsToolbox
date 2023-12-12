function [theta] = ANGLE_BETWEEN(V1, V2)
% [theta] = ANGLE_BETWEEN(V1, V2)
% Find the interior angle between two vectors
    theta = atan2(norm(cross(V1, V2)), dot(V1, V2));
end

