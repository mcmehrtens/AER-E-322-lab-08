% Calculates the moments of inertia and shear centers for lab eight.
function [I_1, I_2, e_1, e_2, e_3, e_4, e_5] = calc_shear_centers
    u = symunit;
    
    % Formula for Moment of Inertia of Specimen 1 & 2
    I = @(h, b, t) (2 * b * t^3 + t * (h - t)^3 + 6 * b * t * h^2) ...
        / 12; % [inch^4]
    
    % Formula for r
    r = @(OD, t) (OD - t) / 2; % [inch]
    
    % Specimen 1 Parameters
    h_1 = 2.43*u.in; % [inch]
    b_1 = 1.456*u.in; % [inch]
    t_1 = 0.08*u.in; % [inch]
    I_1 = I(h_1, b_1, t_1); % [inch^4]
    
    % Specimen 2 Parameters
    h_2 = 0.84*u.in; % [inch]
    b_2 = 0.56*u.in; % [inch]
    t_2 = 0.055*u.in; % [inch]
    I_2 = I(h_2, b_2, t_2); % [inch^4]
    
    % Specimen 3 Parameters
    t_3 = 0.071*u.in; % [inch]
    OD_3 = 1.66*u.in; % [inch]
    theta_0_3 = deg2rad(3.1 / 2); % [rad]
    r_3 = r(OD_3, t_3); % [inch]
    
    % Specimen 4 Parameters
    t_4 = 0.071*u.in; % [inch]
    OD_4 = 1.66*u.in; % [inch]
    theta_0_4 = deg2rad(36.3 / 2); % [rad]
    r_4 = r(OD_4, t_4); % [inch]
    
    % Specimen 5 Parameters
    t_5 = 0.071*u.in; % [inch]
    OD_5 = 1.66*u.in; % [inch]
    theta_0_5 = deg2rad(103.7 / 2); % [rad]
    r_5 = r(OD_5, t_5); % [inch]
    
    % Shear center calculations
    e_c = @(h, b, t, I) -h^2 * b^2 * t / (4 * I); % [inch]
    e_circ = @(r, theta_0) -2 * r ...
        * (cos(theta_0) * (2 * pi - 2 * theta_0) + 2 * sin(theta_0)) ...
        / (2 * pi - 2 * theta_0 + sin(2 * theta_0)); % [inch]
    
    e_1 = e_c(h_1, b_1, t_1, I_1); % [inch]
    e_2 = e_c(h_2, b_2, t_2, I_2); % [inch]
    e_3 = e_circ(r_3, theta_0_3); % [inch]
    e_4 = e_circ(r_4, theta_0_4); % [inch]
    e_5 = e_circ(r_5, theta_0_5); % [inch]
    
    % Print output
    fprintf("I_1 = %g [inch^4] = %g [cm^4]\n" + ...
        "I_2 = %g [inch^4] = %g [cm^4]\n" + ...
        "e_1 = %g [inch] = %g [cm]\n" + ...
        "e_2 = %g [inch] = %g [cm]\n" + ...
        "e_3 = %g [inch] = %g [cm]\n" + ...
        "e_4 = %g [inch] = %g [cm]\n" + ...
        "e_5 = %g [inch] = %g [cm]\n", ...
        double(separateUnits([I_1, unitConvert(I_1, u.cm^4), ...
        I_2, unitConvert(I_2, u.cm^4), e_1, unitConvert(e_1, u.cm), ...
        e_2, unitConvert(e_2, u.cm), e_3, unitConvert(e_3, u.cm), ...
        e_4, unitConvert(e_4, u.cm), e_5, unitConvert(e_5, u.cm)])));
end