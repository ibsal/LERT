function [Radius, Area] = combustionChamberProfile(x, R_combustionchamber, R1, theta1, Rt, T_Throatradius, L, theta2, R_exit, theta_exit, L_nozzle)
    y = zeros(1, length(x));
    for i=1:length(y)
    % Calculate the boundaries for each piece of the function
    boundary1 = -R1 * tan((theta1 * pi / 180) / 2) + (1 / tan(-theta1 * pi / 180)) * (R_combustionchamber + sqrt(Rt^2 - (Rt * sin(theta1 * pi / 180))^2) - T_Throatradius - Rt) + (L - Rt * sin(theta1 * pi / 180));
    
    boundary2 = boundary1 + R1 * sin(theta1 * pi / 180);
    
    boundary3 = L - Rt * sin(theta1 * pi / 180);
    
    boundary4 = L + Rt * sin(theta2 * pi / 180);
    
    exitslope = tan(theta_exit * pi / 180);

    enterslope = tan(theta2 * pi / 180);
    
    XX = [L + Rt * sin(theta2 * pi / 180), L_nozzle + L];
    YY = [-sqrt(Rt^2 - (boundary4 - L)^2) + T_Throatradius + Rt, R_exit];

    nozzleexit = spline(XX, [enterslope YY exitslope]);
    

    % Determine which piece of the function to use based on the value of x
    if x(i) >= 0 && x(i) <= boundary1
        % First piece: y = R_combustionchamber
        y(i) = R_combustionchamber;
    elseif x(i) > boundary1 && x(i) <= boundary2
        % Second piece: y = sqrt(R1^2 - (x - boundary1)^2) + R_combustionchamber - R1
        y(i) = sqrt(R1^2 - (x(i) - boundary1)^2) + R_combustionchamber - R1;
    elseif x(i) > boundary3 && x(i) <= boundary4
        % Third piece: y = -sqrt(Rt^2 - (x - L)^2) + T_Throatradius + Rt
        y(i) = -sqrt(Rt^2 - (x(i) - L)^2) + T_Throatradius + Rt;
    elseif x(i) > boundary2 && x(i) <= boundary3
        % Fourth piece: y = tan(-theta1 * pi / 180) * (x - (L - Rt * sin(theta1 * pi / 180))) - sqrt(Rt^2 - (Rt * sin(theta1 * pi / 180))^2) + T_Throatradius + Rt
        y(i) = tan(-theta1 * pi / 180) * (x(i) - (L - Rt * sin(theta1 * pi / 180))) - sqrt(Rt^2 - (Rt * sin(theta1 * pi / 180))^2) + T_Throatradius + Rt;
    elseif x(i) > boundary4 && x(i) <= XX(2)
        y(i) = ppval(nozzleexit, x(i));
    else
        % If x is outside the defined domain, return NaN (Not a Number)
        y(i) = NaN;
    end
    end

    Radius = y;
    Area = pi .* y.^2;
end