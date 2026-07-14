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

%% Convert ruler position to mirror angle (radians)

thetaForward = atan((forward - xC)/D)/2;
thetaBackward = atan((backward - xC)/D)/2;
thetaMean = (thetaForward + thetaBackward)/2;

%% Linear fits

forwardFit = polyfit(thetaForward, torque, 1);
backwardFit = polyfit(thetaBackward, torque, 1);
meanFit = polyfit(thetaMean, torque, 1);

%% Evaluate fitted lines

thetaPlot = linspace( ...
    min([thetaForward thetaBackward]), ...
    max([thetaForward thetaBackward]), ...
    200);

forwardLine = polyval(forwardFit, thetaPlot);
backwardLine = polyval(backwardFit, thetaPlot);
meanLine = polyval(meanFit, thetaPlot);

%% Plot

figure
hold on

plot(thetaForward, torque, 'bo', 'MarkerFaceColor','b')
plot(thetaBackward, torque, 'ro', 'MarkerFaceColor','r')

plot(thetaPlot, meanLine, 'k-', 'LineWidth',2)

xlabel('\theta (radians)')
ylabel('Torque (g*cm)')
legend('Forward','Backward', ...
       'Mean fit', ...
       'Location','best')

grid on

%% Print calibration constants

fprintf("\nForward:\n");
fprintf("Weight = %.6f * theta + %.6f\n", ...
    forwardFit(1), forwardFit(2));

fprintf("\nBackward:\n");
fprintf("Weight = %.6f * theta + %.6f\n", ...
    backwardFit(1), backwardFit(2));

fprintf("\nMean:\n");
fprintf("Weight = %.6f * theta + %.6f\n", ...
    meanFit(1), meanFit(2));

fprintf("\nSpring calibration coefficients (dyn*cm/rad):\n");
fprintf("Forward : %.6f\n", forwardFit(1));
fprintf("Backward: %.6f\n", backwardFit(1));
fprintf("Mean    : %.6f\n", meanFit(1));