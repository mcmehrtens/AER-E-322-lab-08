% AER E 322 Lab 8 Analysis Script
% Spring 2023 Section 4 Group 2
clear,clc,close all;

% Calculate Angle of Twist and x for Specimen 1
u = symunit;
L_1 = 44.0*u.cm; % [cm]
t_1 = 0.08*u.in; % [in]
t_1 = unitConvert(t_1,u.cm); % [cm]

x_11 = -(11.0*u.cm + t_1 / 2); % [cm]
x_12 = -(14.0*u.cm + t_1 / 2); % [cm]
x_13 = 9.2*u.cm + t_1 / 2; % [cm]
x_14 = 12.8*u.cm + t_1 / 2; % [cm]
x_1 = [x_11 x_12 x_13 x_14]; % [cm]

h_L_1 = [10.4 8.8 17.7 18.7].*u.cm; % [cm]
h_R_1 = [13.1 14.6 4.8 3.4].*u.cm; % [cm]
aot_1 = asind((h_R_1 - h_L_1) ./ L_1)*u.deg; % [deg]
err_1 = asind((9.0*u.cm - 14.0*u.cm) / L_1)*u.deg; % [deg]
aot_1 = aot_1 - err_1; % [deg]

% Calculate Angle of Twist and x for Specimen 2
L_2 = 30.0*u.cm; % [cm]
t_2 = 0.055*u.in; % [in]
t_2 = unitConvert(t_2,u.cm); % [cm]

x_21 = -(3.1*u.cm + t_2 / 2); % [cm]
x_22 = -(9.0*u.cm + t_2 / 2); % [cm]
x_23 = 3.3*u.cm + t_2 / 2; % [cm]
x_24 = 7.5*u.cm + t_2 / 2; % [cm]
x_2 = [x_21 x_22 x_23 x_24]; % [cm]

h_L_2 = [9.9 8.1 11.6 12.6].*u.cm; % [cm]
h_R_2 = [11.4 13.0 9.5 8.3].*u.cm; % [cm]
aot_2 = asind((h_R_2 - h_L_2) ./ L_2)*u.deg; % [deg]

% Calculate Angle of Twist and x for Specimen 3
L_3 = 30.0*u.cm; % [cm]
t_3 = 0.071*u.in; % [in]
t_3 = unitConvert(t_3,u.cm); % [cm]
OD_3 = 1.66*u.in; % [in]
OD_3 = unitConvert(OD_3,u.cm); % [cm]

x_31 = -(4.9*u.cm + OD_3 / 2); % [cm]
x_32 = -(7.4*u.cm + OD_3 / 2); % [cm]
x_33 = 9.8*u.cm + t_3 - OD_3 / 2; % [cm]
x_34 = 5.6*u.cm + t_3 - OD_3 / 2; % [cm]
x_3 = [x_31 x_32 x_33 x_34]; % [cm]

h_L_3 = [7.7 6.9 12.1 11.5].*u.cm; % [cm]
h_R_3 = [12.5 13.4 8.2 9.1].*u.cm; % [cm]
aot_3 = asind((h_R_3 - h_L_3) ./ L_3)*u.deg; % [deg]

% Calculate Angle of Twist and x for Specimen 4
L_4 = 30.0*u.cm; % [cm]
t_4 = 0.071*u.in; % [in]
t_4 = unitConvert(t_4,u.cm); % [cm]
OD_4 = 1.66*u.in; % [in]
OD_4 = unitConvert(OD_4,u.cm); % [cm]

x_41 = -(5.0*u.cm + OD_4 / 2); % [cm]
x_42 = -(7.6*u.cm + OD_4 / 2); % [cm]
x_43 = 7.5*u.cm + t_4 - OD_4 / 2; % [cm]
x_44 = 4.5*u.cm + t_4 - OD_4 / 2; % [cm]
x_4 = [x_41 x_42 x_43 x_44]; % [cm]

h_L_4 = [9.8 8.9 15.3 14.2].*u.cm; % [cm]
h_R_4 = [13.4 14.7 6.1 7.7].*u.cm; % [cm]
aot_4 = asind((h_R_4 - h_L_4) ./ L_4)*u.deg; % [deg]

% Calculate Angle of Twist and x for Specimen 5
L_5 = 30.5*u.cm; % [cm]
t_5 = 0.071*u.in; % [in]
t_5 = unitConvert(t_5,u.cm); % [cm]
OD_5 = 1.66*u.in; % [in]
OD_5 = unitConvert(OD_5,u.cm); % [cm]

x_51 = -(4.5*u.cm + OD_5 / 2); % [cm]
x_52 = -(10.8*u.cm + OD_5 / 2); % [cm]
x_53 = 5.6*u.cm + t_5 - OD_5 / 2; % [cm]
x_54 = 10.6*u.cm + t_5 - OD_5 / 2; % [cm]
x_5 = [x_51 x_52 x_53 x_54]; % [cm]

h_L_5 = [10.0 7.6 15.0 17.2]*u.cm; % [cm]
h_R_5 = [13.8 16.4 7.2 5.8]*u.cm; % [cm]
aot_5 = asind((h_R_5 - h_L_5) / L_5)*u.deg; % [deg]

fprintf("Measured Shear Centers:\n")

% Plot and Calculate Shear Center of Specimen 1
c_1 = polyfit(double(separateUnits(x_1)), double(separateUnits(aot_1)), 1);
lobf_x_1 = linspace(double(separateUnits(x_1(2))), ...
    double(separateUnits(x_1(end))), ...
    1000); % [cm]
lobf_y_1 = c_1(1) * lobf_x_1 + c_1(2); % [deg]
e_1 = (-c_1(2) / c_1(1))*u.cm; % [cm]

fprintf("e_1 = %g [cm]\n", separateUnits(e_1));

figure(1);
scatter(separateUnits(x_1), separateUnits(aot_1), ...
    "DisplayName", "Measurements");
hold on;
plot(lobf_x_1, lobf_y_1, "DisplayName", "Line of Best Fit");
scatter(separateUnits(e_1), 0, "filled", "DisplayName", "Shear Center");
hold off;
title("Angle of Twist vs. Distance from Reference Center");
xlabel("x (cm)");
ylabel("theta (deg)");
legend;
grid on;

% Plot and Calculate Shear Center of Specimen 2
c_2 = polyfit(double(separateUnits(x_2)), double(separateUnits(aot_2)), 1);
lobf_x_2 = linspace(double(separateUnits(x_2(2))), ...
    double(separateUnits(x_2(end))), ...
    1000); % [cm]
lobf_y_2 = c_2(1) * lobf_x_2 + c_2(2); % [deg]
e_2 = (-c_2(2) / c_2(1))*u.cm; % [cm]

fprintf("e_2 = %g [cm]\n", separateUnits(e_2));

figure(2);
scatter(separateUnits(x_2), separateUnits(aot_2), ...
    "DisplayName", "Measurements");
hold on;
plot(lobf_x_2, lobf_y_2, "DisplayName", "Line of Best Fit");
scatter(separateUnits(e_2), 0, "filled", "DisplayName", "Shear Center");
hold off;
title("Angle of Twist vs. Distance from Reference Center");
xlabel("x (cm)");
ylabel("theta (deg)");
legend;
grid on;

% Plot and Calculate Shear Center of Specimen 3
c_3 = polyfit(double(separateUnits(x_3)), double(separateUnits(aot_3)), 1);
lobf_x_3 = linspace(double(separateUnits(x_3(2))), ...
    double(separateUnits(x_3(end - 1))), ...
    1000); % [cm]
lobf_y_3 = c_3(1) * lobf_x_3 + c_3(2); % [deg]
e_3 = (-c_3(2) / c_3(1))*u.cm; % [cm]

fprintf("e_3 = %g [cm]\n", separateUnits(e_3));

figure(3);
scatter(separateUnits(x_3), separateUnits(aot_3), ...
    "DisplayName", "Measurements");
hold on;
plot(lobf_x_3, lobf_y_3, "DisplayName", "Line of Best Fit");
scatter(separateUnits(e_3), 0, "filled", "DisplayName", "Shear Center");
hold off;
title("Angle of Twist vs. Distance from Reference Center");
xlabel("x (cm)");
ylabel("theta (deg)");
legend;
grid on;

% Plot and Calculate Shear Center of Specimen 4
c_4 = polyfit(double(separateUnits(x_4)), double(separateUnits(aot_4)), 1);
lobf_x_4 = linspace(double(separateUnits(x_4(2))), ...
    double(separateUnits(x_4(end - 1))), ...
    1000); % [cm]
lobf_y_4 = c_4(1) * lobf_x_4 + c_4(2); % [deg]
e_4 = (-c_4(2) / c_4(1))*u.cm; % [cm]

fprintf("e_4 = %g [cm]\n", separateUnits(e_4));

figure(4);
scatter(separateUnits(x_4), separateUnits(aot_4), ...
    "DisplayName", "Measurements");
hold on;
plot(lobf_x_4, lobf_y_4, "DisplayName", "Line of Best Fit");
scatter(separateUnits(e_4), 0, "filled", "DisplayName", "Shear Center");
hold off;
title("Angle of Twist vs. Distance from Reference Center");
xlabel("x (cm)");
ylabel("theta (deg)");
legend;
grid on;

% Plot and Calculate Shear Center of Specimen 5
c_5 = polyfit(double(separateUnits(x_5)), double(separateUnits(aot_5)), 1);
lobf_x_5 = linspace(double(separateUnits(x_5(2))), ...
    double(separateUnits(x_5(end))), ...
    1000); % [cm]
lobf_y_5 = c_5(1) * lobf_x_5 + c_5(2); % [deg]
e_5 = (-c_5(2) / c_5(1))*u.cm; % [cm]

fprintf("e_5 = %g [cm]\n", separateUnits(e_5));

figure(5);
scatter(separateUnits(x_5), separateUnits(aot_5), ...
    "DisplayName", "Measurements");
hold on;
plot(lobf_x_5, lobf_y_5, "DisplayName", "Line of Best Fit");
scatter(separateUnits(e_5), 0, "filled", "DisplayName", "Shear Center");
hold off;
title("Angle of Twist vs. Distance from Reference Center");
xlabel("x (cm)");
ylabel("theta (deg)");
legend;
grid on;

% Calculate moments of inertia and theoretical shear centers
fprintf("\nTheoretical Shear Centers:\n")
[I_1, I_2, e_1_th, e_2_th, e_3_th, e_4_th, e_5_th] = calc_shear_centers;
I_1 = unitConvert(I_1, u.m^4); % [m^4]
I_2 = unitConvert(I_2, u.m^4); % [m^4]
e_1_th = unitConvert(e_1_th, u.cm); % [cm]
e_2_th = unitConvert(e_2_th, u.cm); % [cm]
e_3_th = unitConvert(e_3_th, u.cm); % [cm]
e_4_th = unitConvert(e_4_th, u.cm); % [cm]
e_5_th = unitConvert(e_5_th, u.cm); % [cm]

% Calculate relative errors
err = @(meas, theor) abs((theor - meas) / theor) * 100; % [%]
e_1_err = err(e_1, e_1_th);
e_2_err = err(e_2, e_2_th);
e_3_err = err(e_3, e_3_th);
e_4_err = err(e_4, e_4_th);
e_5_err = err(e_5, e_5_th);

fprintf("\nRelative Errors:\n" + ...
    "e_1_err = %g%%\n" + ...
    "e_2_err = %g%%\n" + ...
    "e_3_err = %g%%\n" + ...
    "e_4_err = %g%%\n" + ...
    "e_5_err = %g%%\n", ...
    [e_1_err, e_2_err, e_3_err, e_4_err, e_5_err]);

% General Equations of Shear Flow
tau_12 = @(P, h, I, s_1) P * h * s_1 / (2 * I); % [Pa]
tau_24 = @(P, h, I, b, s_2) tau_12(P, h, I, b) ...
    - P / (2 * I) * (s_2^2 - h * s_2); % [Pa]
tau_45 = @(P, h, I, b, s_4) tau_24(P, h, I, b, h) ...
    - P * h * s_4 / (2 * I); % [Pa]

% Calculate shear flow in specimen 1
fprintf("\nShear Flow Calculations for Specimen 1:\n")
P_1 = unitConvert(100*u.g, u.kg) * (9.81*u.m/u.s^2); % [N]
h_1 = unitConvert(2.43*u.in, u.m); % [m]
b_1 = unitConvert(1.456*u.in, u.m); % [m]

syms s;
tau_12_1 = unitConvert(tau_12(P_1, h_1, I_1, s), u.kPa); % [kPa]
tau_24_1 = unitConvert(tau_24(P_1, h_1, I_1, b_1, s), u.kPa); % [kPa]
tau_45_1 = unitConvert(tau_45(P_1, h_1, I_1, b_1, s), u.kPa); % [kPa]
tau_1_1 = unitConvert(subs(tau_12_1, s, 0), u.kPa); % [kPa]
tau_2_1 = unitConvert(subs(tau_12_1, s, b_1), u.kPa); % [kPa]
tau_3_1 = unitConvert(subs(tau_24_1, s, h_1 / 2), u.kPa); % [kPa]
tau_4_1 = unitConvert(subs(tau_24_1, s, h_1), u.kPa); % [kPa]
tau_5_1 = unitConvert(subs(tau_45_1, s, b_1), u.kPa); % [kPa]

fprintf("tau_12_1 = "); disp(vpa(tau_12_1, 4));
fprintf("tau_24_1 = "); disp(vpa(tau_24_1, 4));
fprintf("tau_45_1 = "); disp(vpa(tau_45_1, 4));
fprintf("tau_1_1 = "); disp(vpa(tau_1_1, 4));
fprintf("tau_2_1 = "); disp(vpa(tau_2_1, 4));
fprintf("tau_3_1 = "); disp(vpa(tau_3_1, 4));
fprintf("tau_4_1 = "); disp(vpa(tau_4_1, 4));
fprintf("tau_5_1 = "); disp(vpa(tau_5_1, 4));

% Calculate shear flow in specimen 2
fprintf("Shear Flow Calculations for Specimen 2:\n")
P_2 = unitConvert(200*u.g, u.kg) * (9.81*u.m/u.s^2); % [N]
h_2 = unitConvert(0.84*u.in, u.m); % [m]
b_2 = unitConvert(0.56*u.in, u.m); % [m]

syms s;
tau_12_2 = unitConvert(tau_12(P_2, h_2, I_2, s), u.kPa); % [kPa]
tau_24_2 = unitConvert(tau_24(P_2, h_2, I_2, b_2, s), u.kPa); % [kPa]
tau_45_2 = unitConvert(tau_45(P_2, h_2, I_2, b_2, s), u.kPa); % [kPa]
tau_1_2 = unitConvert(subs(tau_12_2, s, 0), u.kPa); % [kPa]
tau_2_2 = unitConvert(subs(tau_12_2, s, b_2), u.kPa); % [kPa]
tau_3_2 = unitConvert(subs(tau_24_2, s, h_2 / 2), u.kPa); % [kPa]
tau_4_2 = unitConvert(subs(tau_24_2, s, h_2), u.kPa); % [kPa]
tau_5_2 = unitConvert(subs(tau_45_2, s, b_2), u.kPa); % [kPa]

fprintf("tau_12_2 = "); disp(vpa(tau_12_2, 4));
fprintf("tau_24_2 = "); disp(vpa(tau_24_2, 4));
fprintf("tau_45_2 = "); disp(vpa(tau_45_2, 4));
fprintf("tau_1_2 = "); disp(vpa(tau_1_2, 4));
fprintf("tau_2_2 = "); disp(vpa(tau_2_2, 4));
fprintf("tau_3_2 = "); disp(vpa(tau_3_2, 4));
fprintf("tau_4_2 = "); disp(vpa(tau_4_2, 4));
fprintf("tau_5_2 = "); disp(vpa(tau_5_2, 4));