function rulerInches = projectToRulerAxis(pixelPosition,rulerCalibration)

x = pixelPosition(:,1);
y = pixelPosition(:,2);

% Warn if the laser drifts vertically from the calibrated camera axis.
verticalError = abs(y-rulerCalibration.cameraY);

if verticalError > 20
    warning("Laser is %.1f pixels away from the calibrated camera axis.",verticalError)
end

% Convert horizontal coordinate to ruler inches.
rulerInches = ppval(rulerCalibration.pixelToInchSpline,x);

end
