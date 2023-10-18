function [B] = OA2BA(O)
    % Solve Barker's Equation (cubic) using Cardano's Method
    % (Vallado 1027)
    P = 0;
    Q = 3;
    R = -3*O;
    a = (1/3)*(3*Q - P^2); % = 3
    b = (1/27)*(2*P^3 - 9*P*Q + 27*R); % = R
    del = a^3/27 + b^2/4; % always positive

    % Calculate the one real root
    x = (-b/2 + sqrt(del))^(1/3) + (-b/2 - sqrt(del))^(1/3);

    B = x;
end