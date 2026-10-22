clc
close all
clear

% Read shared batch paths; the fitted axis and spline method below remain
% unchanged. The same saved result is reused for all batch trials.
calibrationDir = fileparts(mfilename('fullpath'));
addpath(fileparts(calibrationDir));
config = loadBatchConfig();

%% SETUP
rulerValues = [8 10 12 14 17.5 22 24 26];
numPoints = numel(rulerValues);

originValue = 17.5;          % Physical value chosen as xC

% Refresh the displayed reference whenever calibration is explicitly started.
% Use repeat 1 of the first configured group, just as save_frames.m does.
% Copy the extracted JPEG directly so its pixels match the tracking images.
referenceFrame = fullfile(config.paths.frames, config.trials(1) + "_frames", '001.jpg');
assert(isfile(referenceFrame), ...
    'Missing first batch frame: %s. Run save_frames.m before calibration.', referenceFrame);
referenceDir = fileparts(config.paths.calibration_reference);
if ~isfolder(referenceDir), mkdir(referenceDir); end
copyfile(referenceFrame, config.paths.calibration_reference, 'f');
imageData = imread(config.paths.calibration_reference);

%% COLLECT CLICKS

calibrationFigure = figure('Name','Ruler calibration','NumberTitle','off');
calibrationAxes = axes('Parent',calibrationFigure);
imshow(imageData,'Parent',calibrationAxes)
% Leave room for the instructions above and the confirmation button below.
set(calibrationAxes,'Position',[0.05 0.14 0.9 0.70]);
pointHandles = cell(numPoints,1);
for pointIndex = 1:numPoints
    title(calibrationAxes, {
        sprintf('Click the %g-inch ruler mark (%d of %d)', ...
                rulerValues(pointIndex),pointIndex,numPoints), ...
        'Each point stays labeled; drag markers to correct their positions.'});
    % Image Processing Toolbox point ROIs keep previous clicks visible and
    % editable. Labels identify the physical reading, not the click number.
    pointHandles{pointIndex} = drawpoint(calibrationAxes, ...
        'Color','yellow','Label',sprintf('%g in',rulerValues(pointIndex)), ...
        'Deletable',false);
    if ~isgraphics(calibrationFigure) || isempty(pointHandles{pointIndex}) ...
            || ~isvalid(pointHandles{pointIndex})
        error('Calibration cancelled; the saved calibration was not changed.');
    end
end

% During review, a clicked point becomes the keyboard target (shown in cyan).
% Keep listener handles alive while uiwait runs; dragging remains available.
calibrationFigure.UserData = struct('points',{pointHandles}, ...
    'activePoint',pointHandles{1},'imageSize',size(imageData));
pointListeners = cell(numPoints,1);
for pointIndex = 1:numPoints
    pointListeners{pointIndex} = addlistener(pointHandles{pointIndex}, ...
        'ROIClicked',@(source,~) selectCalibrationPoint(calibrationFigure,source));
end
selectCalibrationPoint(calibrationFigure,pointHandles{1});
calibrationFigure.WindowKeyPressFcn = @nudgeCalibrationPoint;

% Read positions only after review, so both dragging and nudging update the fit.
title(calibrationAxes, {'Review: click a point to select it (cyan), or drag it.', ...
    'Arrow keys: 1 pixel. Shift + arrow: 10 pixels.', ...
    'Click Accept points when ready. Close the window to cancel.'});
uicontrol(calibrationFigure,'Style','pushbutton','String','Accept points', ...
    'Units','normalized','Position',[0.39 0.03 0.22 0.06], ...
    'Callback',@(~,~) uiresume(calibrationFigure));
uiwait(calibrationFigure);
if ~isgraphics(calibrationFigure)
    error('Calibration cancelled; the saved calibration was not changed.');
end
clickedPoints = zeros(numPoints,2);
for pointIndex = 1:numPoints
    clickedPoints(pointIndex,:) = pointHandles{pointIndex}.Position;
end
x = clickedPoints(:,1);
y = clickedPoints(:,2);

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

% Honor a configured calibration destination, including a new batch folder.
outputDir = fileparts(config.paths.ruler_calibration);
if ~isfolder(outputDir), mkdir(outputDir); end
save(config.paths.ruler_calibration,"rulerCalibration")

function selectCalibrationPoint(fig, point)
% Keep one explicit keyboard target, independent of which ROI has GUI focus.
    state = fig.UserData;
    for index = 1:numel(state.points)
        state.points{index}.Color = 'yellow';
    end
    point.Color = 'cyan';
    state.activePoint = point;
    fig.UserData = state;
end

function nudgeCalibrationPoint(fig, event)
% Image y increases downward: up subtracts one pixel from the y coordinate.
    switch event.Key
        case 'leftarrow',  direction = [-1 0];
        case 'rightarrow', direction = [1 0];
        case 'uparrow',    direction = [0 -1];
        case 'downarrow',  direction = [0 1];
        otherwise, return;
    end
    step = 1;
    if any(strcmp(event.Modifier,'shift')), step = 10; end
    state = fig.UserData;
    point = state.activePoint;
    position = point.Position + step*direction;
    % Prevent nudging a marker beyond the image's pixel centers.
    point.Position = max([1 1], min(position, state.imageSize([2 1])));
end
