clc
clear
close all

%% Calibration constants

D = 48.5;      % Distance from camera to mirror (in)
xC = 18;       % Closest ruler position to the camera (in)
g = 980;       % Acceleration due to gravity (cm/s^2)
a = 5.95/2;      % Distance from mirror to spring (cm)

hook_mass = 0.148;      % g
added_mass = 0.182;     % g

% Total applied weights
mass = hook_mass + (0:8)*added_mass;

% Torque
torque = mass * g * 2 * a;

%% Measured ruler positions (inches)

forward = [5, 8.375, 11.375, 15, 18, 21.875, 24.75, 28.125, 32];
backward = [6.625, 9.75, 11.625, 15.375, 18.25, 21.75, 24.75, 28.375, 31.875];

% Join both
position = [forward backward];
torque = [torque torque];

%% Convert ruler position to mirror angle (radians)
theta = atan((position - xC)/D)/2;

%% Linear fits
fit = polyfit(theta, torque, 1);

thetaPlot = linspace( ...
    min([theta ]), ...
    max([theta ]), ...
    200);

line = polyval(fit, thetaPlot);

%% Plot

figure
hold on

plot(theta, torque, 'bo', 'MarkerFaceColor','b')

plot(thetaPlot, line, 'k-', 'LineWidth',2)

xlabel('\theta (radians)')
ylabel('Torque (g*cm)')
legend('theta', ...
       'fit', ...
       'Location','best')

grid on

%% Print calibration constants

fprintf("Weight = %.6f * theta + %.6f\n", ...
    fit(1), fit(2));

fprintf("\nSpring calibration coefficients (dyn*cm/rad):\n");
fprintf("%.6f\n", fit(1));