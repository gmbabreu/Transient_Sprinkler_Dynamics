function limits = calibrateLaserMask(configFile, imageFile, center)
%CALIBRATELASERMASK Click the red laser, preview its mask, then save JSON limits.
% Run from the editor, or add this folder to the path and call calibrateLaserMask.
% Optional configFile/imageFile select a different configuration/reference image.
% Supplying center=[x y] returns a suggestion without opening a GUI or saving.
    if nargin < 1 || isempty(configFile)
        configFile = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'experiment_config.json');
    end
    addpath(fileparts(fileparts(mfilename('fullpath'))));
    [config, trialPath] = loadBatchConfig(configFile);
    if nargin < 2 || isempty(imageFile)
        imageFile = config.paths.calibration_reference;
        if isfield(config, 'calibration') && isfield(config.calibration, 'reference_trial') ...
                && strlength(string(config.calibration.reference_trial)) > 0
            imageFile = fullfile(trialPath(config.paths.frames, ...
                string(config.calibration.reference_trial), "_frames"), '001.jpg');
        end
    end
    assert(isfile(imageFile), 'Missing calibration image: %s', imageFile);
    img = imread(imageFile);
    radius = config.video.circle_radius_pixels;
    if nargin >= 3
        limits = suggestLimits(img, center, radius);
        return;
    end
    originalText = fileread(configFile);
    limits = [];
    choice = questdlg('Preview the saved mask limits or select a new laser sample?', ...
        'Laser mask calibration', 'Preview current', 'Calibrate', 'Cancel', 'Preview current');
    if isempty(choice) || strcmp(choice, 'Cancel'), return; end
    fig = figure('Name', 'Laser mask calibration', 'NumberTitle', 'off', 'WindowState', 'maximized');
    if strcmp(choice, 'Preview current')
        showMask(fig, img, config.video.laser_ycbcr_limits, 'Current mask');
        fprintf('Current Y/Cb/Cr limits: %s\n', jsonencode(config.video.laser_ycbcr_limits));
        choice = questdlg('Current settings are unchanged. Continue to recalibrate?', ...
            'Current mask preview', 'Calibrate', 'Close', 'Close');
        if ~strcmp(choice, 'Calibrate') || ~isgraphics(fig)
            if isgraphics(fig), close(fig); end
            return;
        end
    end
    while isgraphics(fig)
        clf(fig);
        ax = axes('Parent', fig);
        imshow(img, 'Parent', ax);
        title(ax, {'Click near the red laser dot to zoom in.', ...
            'Close the window to cancel.'});
        point = drawpoint(ax, 'Color', 'cyan');
        if ~isgraphics(fig) || isempty(point), return; end
        center = point.Position;
        delete(point);
        halfWidth = max(50, 8*radius);
        xlim(ax, [max(.5, center(1)-halfWidth), min(size(img,2)+.5, center(1)+halfWidth)]);
        ylim(ax, [max(.5, center(2)-halfWidth), min(size(img,1)+.5, center(2)+halfWidth)]);
        title(ax, 'Drag/resize the circle around the laser, then double-click it.');
        region = drawcircle(ax, 'Center', center, 'Radius', radius, 'Color', 'cyan');
        wait(region);
        if ~isgraphics(fig) || ~isvalid(region), return; end
        try
            proposed = suggestLimits(img, region.Center, region.Radius);
        catch err
            errordlg(err.message, 'Choose a clearer laser sample', 'modal');
            continue;
        end
        showMask(fig, img, proposed, 'Proposed mask');
        fprintf('Suggested Y/Cb/Cr limits: %s\n', jsonencode(proposed));
        choice = questdlg('Inspect the whole mask for background pixels before saving.', ...
            'Laser mask calibration', 'Save', 'Pick again', 'Cancel', 'Cancel');
        if strcmp(choice, 'Pick again'), continue; end
        if ~strcmp(choice, 'Save') || ~isgraphics(fig)
            if isgraphics(fig), close(fig); end
            return;
        end
        % Replace only this value: decoding/re-encoding the entire config would
        % change numeric-leading group names such as 1000r into MATLAB fields.
        assert(strcmp(fileread(configFile), originalText), ...
            'Configuration changed during calibration; rerun to avoid losing edits.');
        pattern = '"laser_ycbcr_limits"\s*:\s*\[\s*\[[^\]]*\]\s*,\s*\[[^\]]*\]\s*,\s*\[[^\]]*\]\s*\]';
        assert(numel(regexp(originalText, pattern)) == 1, 'Expected one laser_ycbcr_limits matrix.');
        updated = regexprep(originalText, pattern, ['"laser_ycbcr_limits": ' jsonencode(proposed)]);
        parent = char(fileparts(configFile));
        if isempty(parent), parent = pwd; end
        temporary = string(tempname(parent)) + ".json";
        fid = fopen(temporary, 'w');
        assert(fid ~= -1, 'Cannot write calibration settings.');
        cleanup = onCleanup(@() fclose(fid));
        bytes = unicode2native(updated, 'UTF-8');
        written = fwrite(fid, bytes, 'uint8');
        assert(written == numel(bytes), 'Incomplete write; original configuration was not replaced.');
        clear cleanup;
        [ok, message] = movefile(temporary, configFile, 'f');
        assert(ok, '%s', message);
        limits = proposed;
        close(fig);
        fprintf('Saved laser limits to %s. Existing trial CSVs were not regenerated.\n', configFile);
        return;
    end
end

function showMask(fig, img, limits, label)
    [BW, masked] = ycbcrMask(img, limits);
    clf(fig);
    ax = subplot(1,2,1, 'Parent', fig);
    imshow(img, 'Parent', ax); title(ax, 'Reference image');
    ax = subplot(1,2,2, 'Parent', fig);
    imshow(masked, 'Parent', ax);
    title(ax, sprintf('%s: %d pixels. Only the laser should remain.', label, nnz(BW)));
end

function limits = suggestLimits(img, center, radius)
% Separate the red spot from its local background using Cr (red chroma).
% Leave Y/Cb unrestricted so moderate brightness changes do not reject the dot.
    validateattributes(center, {'numeric'}, {'numel',2,'real','finite'});
    validateattributes(radius, {'numeric'}, {'scalar','real','finite','positive'});
    [x,y] = meshgrid(1:size(img,2), 1:size(img,1));
    distance = hypot(x-center(1), y-center(2));
    colors = rgb2ycbcr(img);
    cr = double(colors(:,:,3));
    spot = sort(cr(distance <= radius));
    background = sort(cr(distance >= 2*radius & distance <= 4*radius));
    assert(numel(spot) >= 5 && numel(background) >= 5, 'Select the dot away from the image edge.');
    spotRed = spot(ceil(.9*numel(spot)));
    backgroundRed = background(ceil(.99*numel(background)));
    assert(spotRed > backgroundRed + 5, ...
        'Selected spot is not clearly redder than its surroundings. Recenter or resize the circle.');
    cutoff = ceil((spotRed + backgroundRed)/2);
    limits = [0 255; 0 255; cutoff 255];
end
