function rulerInches = projectToRulerAxis(pixelPosition,rulerCalibration)

x = pixelPosition(:,1);
y = pixelPosition(:,2);

m = rulerCalibration.axisSlope;
cameraY = rulerCalibration.cameraY;

% Move the detected laser horizontally onto the camera axis.
xProjected = x + (cameraY - y)/m;

% Convert projected x-coordinate to ruler inches.
rulerInches = ppval(rulerCalibration.pixelToInchSpline,xProjected);

end
