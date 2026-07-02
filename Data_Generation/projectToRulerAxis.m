function rulerInches = projectToRulerAxis(pixelPosition, rulerCalibration)
%PROJECTTORULERAXIS Project image coordinates onto the fitted ruler axis and convert to inches.

if size(pixelPosition, 2) ~= 2
    error("pixelPosition must be an N-by-2 array of [x y] coordinates.")
end

origin = rulerCalibration.axisOrigin;
axisDirection = rulerCalibration.axisDirection;
axisDirection = axisDirection(:)';

projectedDistance = (pixelPosition - origin) * axisDirection';

if isfield(rulerCalibration, "pixelToInchCoefficients")
    rulerInches = polyval(rulerCalibration.pixelToInchCoefficients, projectedDistance);
else
    rulerInches = projectedDistance;
end

end
