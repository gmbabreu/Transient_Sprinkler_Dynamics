clc
close all
clear

% Calibration settings.
numPoints = 5;
rulerValues = [34 33 32 31 30];

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

% Click the ruler points in image space.
figure("Name", "Ruler Calibration", "Color", "w");
imshow(imageData);
title("Click " + numPoints + " points along the ruler, in order from one end to the other.");

[clickedX, clickedY] = ginput(numPoints);
clickedPoints = [clickedX, clickedY];

if numel(clickedX) < 2
    error("At least two points are required to fit the ruler axis.")
end

axisOrigin = mean(clickedPoints, 1);
centeredPoints = clickedPoints - axisOrigin;
[~, ~, singularVectors] = svd(centeredPoints, 'econ');
axisDirection = singularVectors(:, 1)';

if dot(clickedPoints(end, :) - clickedPoints(1, :), axisDirection) < 0
    axisDirection = -axisDirection;
end

% Project clicks onto the fitted ruler axis and fit pixels to inches.
projectedDistances = (clickedPoints - clickedPoints(1, :)) * axisDirection';
pixelToInchCoefficients = polyfit(projectedDistances, rulerValues, 1);

% Visualize the fitted axis on the image.
hold on
lineHalfLength = max(abs(projectedDistances)) + 50;
lineStart = axisOrigin - lineHalfLength * axisDirection;
lineEnd = axisOrigin + lineHalfLength * axisDirection;
plot(clickedPoints(:, 1), clickedPoints(:, 2), 'yo', 'MarkerSize', 8, 'LineWidth', 1.5)
plot([lineStart(1), lineEnd(1)], [lineStart(2), lineEnd(2)], 'c-', 'LineWidth', 2)
for pointIndex = 1:size(clickedPoints, 1)
    text(clickedPoints(pointIndex, 1) + 5, clickedPoints(pointIndex, 2), sprintf('%.0f', rulerValues(pointIndex)), ...
        'Color', 'w', 'FontWeight', 'bold')
end
hold off

% Plot the calibration fit in pixels vs inches.
figure("Name", "Projected Ruler Distances", "Color", "w");
plot(projectedDistances, rulerValues, 'ko', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
hold on
fitDistances = linspace(min(projectedDistances), max(projectedDistances), 100);
plot(fitDistances, polyval(pixelToInchCoefficients, fitDistances), 'c-', 'LineWidth', 2)
hold off
xlabel("Distance along ruler axis (pixels)")
ylabel("Ruler value (inches)")
grid on

rulerCalibration = struct();
rulerCalibration.imagePath = imagePath;
rulerCalibration.clickedPoints = clickedPoints;
rulerCalibration.axisOrigin = axisOrigin;
rulerCalibration.axisDirection = axisDirection;
rulerCalibration.projectedDistances = projectedDistances;
rulerCalibration.rulerValues = rulerValues;
rulerCalibration.pixelToInchCoefficients = pixelToInchCoefficients;

% Print the final pixel-to-inch relationship.
pixelToInchSlope = pixelToInchCoefficients(1);
pixelToInchIntercept = pixelToInchCoefficients(2);

fprintf("\nCalibration fit:\n");
fprintf("  inches = %.8f * pixelDistance + %.8f\n", pixelToInchSlope, pixelToInchIntercept);
fprintf("  effective inches per pixel along ruler axis: %.8f\n", pixelToInchSlope);

% Save the calibration next to the frame folder.
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
