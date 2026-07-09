clc
close all
clear

%% SETUP 

% Calibration settings.
numPoints = 5;
rulerValues = [10, 14, 18, 22, 26]; % Ruler values in inches corresponding to clicked points

if numel(rulerValues) ~= numPoints
    error("rulerValues must have the same number of entries as numPoints.")
end

% Load a frame image to calibrate against.
if isfile("1500f1_frames/001.jpg")
    imagePath = "1500f1_frames/001.jpg";
else
    [fileName, folderName] = uigetfile({"*.jpg;*.jpeg;*.png;*.tif;*.tiff;*.bmp", "Image files"}, ...
        "Select a frame image for ruler calibration");
    if isequal(fileName, 0)
        error("No calibration image selected.")
    end
    imagePath = fullfile(folderName, fileName);
end

imageData = imread(imagePath);

% Optionally undistort the calibration image if camera calibration data is available.
cameraCalFile = "cameraParams.mat";
useUndistort = isfile(cameraCalFile);
if useUndistort
    % Rectify lens distortion before the user clicks ruler points.
    cameraData = load(cameraCalFile, "cameraParams");
    cameraParams = cameraData.cameraParams;
    imageData = undistortImage(imageData, cameraParams);
end

%% COLLECT CLICKS

% Click the ruler points in image space.
figure("Name", "Ruler Calibration", "Color", "w");
imshow(imageData);
title("Click " + numPoints + " points along the ruler, in order from one end to the other.");

[clickedX, clickedY] = ginput(numPoints);
clickedPoints = [clickedX, clickedY];

if numel(clickedX) < 2
    error("At least two points are required to fit the ruler axis.")
end


%% FIT RULER AXIS

% Fit the ruler axis with polyfit
% This gives a simple line model y = mx + b in image coordinates.
lineCoefficients = polyfit(clickedPoints(:, 1), clickedPoints(:, 2), 1);
axisSlope = lineCoefficients(1);
axisIntercept = lineCoefficients(2);

%% CAMERA AXIS

title("Click the origin position")

[~,cameraY] = ginput(1);
cameraY = cameraY(1);

%% HORIZONTAL PROJECTION CALIBRATION

% Project every ruler mark horizontally onto the camera axis.
% Fit a smooth monotone spline from horizontal pixel coordinate to ruler inches.
pixelToInchSpline = pchip(clickedX,rulerValues);

predictedInches = ppval(pixelToInchSpline,clickedX);

calibrationRmse = sqrt(mean((predictedInches(:)-rulerValues(:)).^2));

%% VISUALIZE AND SAVE

% Visualize the fitted axis on the image.
figure("Name","Calibration","Color","w");
imshow(imageData)
hold on

% Draw the fitted ruler line.
xPlot = linspace(1,size(imageData,2),400);
plot(xPlot,...
     axisSlope*xPlot + axisIntercept,...
     'c-','LineWidth',2)

% Draw the horizontal camera axis.
plot(xPlot,...
     cameraY*ones(size(xPlot)),...
     'r--','LineWidth',2)

% Show the clicked ruler points.
plot(clickedX,clickedY,...
     'yo',...
     'MarkerFaceColor','y',...
     'MarkerSize',8)

% Label the ruler values.
for k = 1:numPoints
    text(clickedX(k)+5,...
         clickedY(k),...
         sprintf("%.0f",rulerValues(k)),...
         "Color","w",...
         "FontWeight","bold");
end

hold off

% Store the calibration outputs that generateData.m needs later.
rulerCalibration = struct();
rulerCalibration.imagePath = imagePath;
rulerCalibration.clickedPoints = clickedPoints;
rulerCalibration.axisSlope = axisSlope;
rulerCalibration.axisIntercept = axisIntercept;
rulerCalibration.rulerValues = rulerValues;
rulerCalibration.pixelToInchSpline = pixelToInchSpline;
rulerCalibration.cameraY = cameraY;

% Print the final pixel-to-inch relationship.
fprintf("\nCalibration fit:\n");
fprintf("  ruler line: y = %.8f*x + %.8f\n",axisSlope,axisIntercept);
fprintf("  pixel-to-inch mapping: PCHIP spline\n");
fprintf("  calibration RMSE = %.6f inches\n",calibrationRmse);

% Save the calibration next to the frame folder.
% The saved MAT file is what generateData.m loads at runtime.
[imageFolder, ~, ~] = fileparts(imagePath);
[parentFolder, frameFolder, ~] = fileparts(imageFolder);
if endsWith(frameFolder, "_frames")
    calibrationName = extractBefore(frameFolder, strlength(frameFolder) - strlength("_frames") + 1) + "_rulerCalibration.mat";
    calibrationPath = fullfile(parentFolder, calibrationName);
else
    calibrationPath = fullfile(imageFolder, "rulerCalibration.mat");
end
save(calibrationPath, "rulerCalibration");

disp("Saved calibration to " + calibrationPath)
