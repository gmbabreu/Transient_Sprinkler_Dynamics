clc
close all
clear

%% SETUP
file = "frame"
rulerValues = [8 10 14 17.5 22 26 28];
numPoints = numel(rulerValues);

originValue = 17.5;          % Physical value chosen as xC

imageData = imread(file + ".jpg");

%% COLLECT CLICKS

figure
imshow(imageData)
title("Click ruler marks in order")

[x,y] = ginput(numPoints);

clickedPoints = [x y];

%% FIT RULER

lineCoefficients = polyfit(x,y,1);

axisSlope = lineCoefficients(1);
axisIntercept = lineCoefficients(2);

axisDirection = [1 axisSlope];
axisDirection = axisDirection/norm(axisDirection);

% Flip so direction agrees with click order
if dot(clickedPoints(end,:)-clickedPoints(1,:),axisDirection)<0
    axisDirection = -axisDirection;
end

%% CHOOSE ORIGIN

originIndex = find(rulerValues==originValue);

if isempty(originIndex)
    error("originValue not found in rulerValues.")
end

% Project this point onto the ruler axis to get the origin point in pixel coordinates
projection = (clickedPoints(originIndex,:) - clickedPoints(1,:))*axisDirection';
originPoint = clickedPoints(1,:) + projection*axisDirection;

%% PROJECT EVERY CLICK ONTO RULER

projectionLength = (clickedPoints-originPoint)*axisDirection';

projectedPoints = originPoint + projectionLength.*axisDirection;

%% TRUE DISTANCES

trueDistances = rulerValues-originValue;

%% FIT NONLINEAR MAPPING

pixelToDistanceSpline = pchip(projectionLength,trueDistances);

predictedDistance = ppval(pixelToDistanceSpline,projectionLength);

%% DISPLAY

figure
imshow(imageData)
hold on

plot(clickedPoints(:,1),clickedPoints(:,2),'yo','MarkerFaceColor','y')

plot(projectedPoints(:,1),projectedPoints(:,2),'r.')

xPlot = linspace(min(x)-50,max(x)+50,300);

plot(xPlot,...
     axisSlope*xPlot+axisIntercept,...
     'c-','LineWidth',2)

hold off

figure

plot(projectionLength,trueDistances,'ko','MarkerFaceColor','k')

hold on

xx = linspace(min(projectionLength),max(projectionLength),300);

plot(xx,ppval(pixelToDistanceSpline,xx),'b','LineWidth',2)

xlabel("Signed pixel distance from origin")
ylabel("True distance from origin (inches)")

grid on

%% SAVE

rulerCalibration.axisSlope = axisSlope;
rulerCalibration.axisIntercept = axisIntercept;
rulerCalibration.axisDirection = axisDirection;

rulerCalibration.originPoint = originPoint;
rulerCalibration.originValue = originValue;

rulerCalibration.pixelToDistanceSpline = pixelToDistanceSpline;

save(file + "_rulerCalibration.mat","rulerCalibration")

