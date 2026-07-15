function inches = projectToRulerAxis(pixelPosition,rulerCalibration)

origin = rulerCalibration.originPoint;
direction = rulerCalibration.axisDirection;

% Signed distance along the ruler
pixelDistance = (pixelPosition-origin)*direction';

% Convert pixel distance to physical distance
distance = ppval(rulerCalibration.pixelToDistanceSpline,pixelDistance);

% Convert to absolute ruler reading
inches = rulerCalibration.originValue + distance;

end