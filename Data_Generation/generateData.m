% Generate individual angle CSVs and timing metadata for one configured batch.
% save_frames.m supplies completed JPEG folders; calibrateRulerAxis.m supplies
% the shared ruler projection. This script also replaces plotData.m's useful
% individual-trial comparison, using audio alignment without averaging.
clc
close all
clear

% Resolve files relative to this script, independently of MATLAB's current folder.
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir, fullfile(scriptDir, 'calibration'));
config = loadBatchConfig();
if ~isfolder(config.paths.data)
    mkdir(config.paths.data);
end

%% Calibration
assert(isfile(config.paths.ruler_calibration), ...
    'Run calibration/calibrateRulerAxis.m first: %s', config.paths.ruler_calibration);
calibrationData = load(config.paths.ruler_calibration, 'rulerCalibration');
% Load one ruler calibration for all listed trials. xC comes from its chosen
% physical ruler origin, while D is the camera/mirror distance in inches.
rulerCalibration = calibrationData.rulerCalibration;
xC = rulerCalibration.originValue;
D = config.video.mirror_distance_inches; % Spring-calibration geometry is separate.

% Clear generated group products before rebuilding them from valid trials.
% Removing the final repeat number gives names such as 1000f. Previous group
% summaries/plots may include a trial that fails this run, so discard them now.
groups = regexprep(config.trials, '[0-9]+$', '');
for group = unique(groups(:))'
    outputDir = fullfile(config.paths.output, group);
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    patterns = {'*_signal.csv', '*_impulses.csv', '*_impulse_diagnostics.png', 'impulse_summary.csv', ...
        'torque.png', 'cumulative_integral.png', 'impulse_comparison.png', 'angular_displacement.png'};
    for p = 1:numel(patterns)
        oldFiles = dir(fullfile(outputDir, patterns{p}));
        for j = 1:numel(oldFiles)
            delete(fullfile(outputDir, oldFiles(j).name));
        end
    end
end
validTrials = struct('id', {}, 'group', {}, 'time', {}, 'angle', {}, 'startup', {});

for trialIndex = 1:numel(config.trials)
    file = config.trials(trialIndex);
    csvPath = fullfile(config.paths.data, file + "_data.csv");
    metadataPath = fullfile(config.paths.data, file + "_metadata.json");
    % Remove old angles before attempting generation; failure must leave no valid CSV.
    if isfile(csvPath)
        delete(csvPath);
    end
    V = [];
    overlayV = [];
    metadata = struct();
    try
        if isfile(metadataPath)
            metadata = jsondecode(fileread(metadataPath));
        end
        % Preserve all saved manual times, fit start and windows on ordinary reruns.
        % Empty arrays mean 'not selected yet'; Python also accepts JSON null.
        % Do not initialize stored integration bounds here: Python derives
        % initial suggestions from the effective event times when needed.
        if ~isfield(metadata, 'manual_event_times_s')
            metadata.manual_event_times_s = struct('startup', [], 'shutdown', []);
        end
        if ~isfield(metadata, 'fit_start_s')
            metadata.fit_start_s = [];
        end
        if ~isfield(metadata, 'integration_windows_s')
            metadata.integration_windows_s = struct('startup', [], 'shutdown', []);
        end
        % Mark incomplete before tracking starts. Only a fully successful run
        % restores valid=true, which Python requires before analyzing the CSV.
        metadata.valid = false;
        metadata.invalid_frames = [];
        metadata.failure_reason = 'Generation incomplete';
        writeMetadata(metadataPath, metadata);

        %% READ, MASK, AND EXTRACT DATA
        folder = fullfile(config.paths.frames, file + "_frames");
        marker = fullfile(folder, '.complete');
        assert(isfile(marker), 'Missing complete extraction; run save_frames.m for %s.', file);
        info = sscanf(fileread(marker), '%f');
        % Reuse extraction's frame rate for timestamps; no second video decode
        % is needed for tracking. Check numbering so a missing JPEG cannot
        % silently truncate the measurement or shift subsequent timestamps.
        assert(numel(info) == 2 && all(isfinite(info)) && info(1) > 0 ...
            && info(1) == fix(info(1)) && info(2) > 0, 'Invalid extraction marker: %s', marker);
        count = info(1);
        fps = info(2);
        expected = compose('%03d.jpg', (1:count)');
        missing = find(~isfile(fullfile(folder, expected)));
        metadata.invalid_frames = missing(:)';
        assert(isempty(missing) && numel(dir(fullfile(folder, '*.jpg'))) == count, ...
            'Incomplete/inconsistent frame extraction; rerun save_frames.m for %s.', file);
        % NaN deliberately survives any missing detection until the whole-trial
        % validity check. Frame 1 corresponds to original video time zero.
        position = NaN(count, 1);
        t = (0:count-1)' / fps;
        x0 = NaN;

        % Save videos.
        if config.video.masked_video
            V = VideoWriter(fullfile(config.paths.data, file + "_maskedData.mp4"), 'MPEG-4');
            V.FrameRate = fps;
            V.Quality = 95;
            open(V);
        end
        if config.video.overlay_video
            overlayV = VideoWriter(fullfile(config.paths.data, file + "_overlayData.mp4"), 'MPEG-4');
            overlayV.FrameRate = fps;
            overlayV.Quality = 95;
            open(overlayV);
        end
        overlayFolder = fullfile(config.paths.data, file + "_overlayFrames");
        if config.video.overlay_frames && ~isfolder(overlayFolder)
            mkdir(overlayFolder);
        end

        for i = 1:count
            % Read the JPEG as an image matrix, then apply the mask defined
            % using the Color Thresholder app (Image Processing Toolbox).
            img = imread(fullfile(folder, expected(i)));
            [BW, maskedImage] = ycbcrMask(img, config.video.laser_ycbcr_limits);
            % Keep the original brightness-weighted centroid of all masked
            % laser pixels, including both x and y for the ruler projection.
            laserImage = sum(maskedImage, 3);
            [laserRows, laserCols] = find(BW);
            laserWeights = double(laserImage(BW));
            if isempty(laserCols) || sum(laserWeights) == 0
                xMid = NaN;
                yMid = NaN;
            else
                xMid = mean(laserCols, 'Weights', laserWeights);
                yMid = mean(laserRows, 'Weights', laserWeights);
            end
            % Do not project missing detections or substitute/interpolate a position.
            if isfinite(xMid) && isfinite(yMid)
                position(i) = projectToRulerAxis([xMid, yMid], rulerCalibration);
            end
            % Retain the first detected ruler position as the angle reference.
            % A missing first frame still invalidates the trial at the end.
            if isnan(x0) && isfinite(position(i))
                x0 = position(i);
            end
            if config.video.masked_video
                writeVideo(V, maskedImage);
            end
            if config.video.overlay_video || config.video.overlay_frames
                % Overlay a hollow black circle at the predicted laser location.
                overlayImage = img;
                if isfinite(xMid) && isfinite(yMid)
                    overlayImage = drawHollowBlackCircle(img, xMid, yMid, ...
                        config.video.circle_radius_pixels, config.video.circle_thickness_pixels);
                end
                if config.video.overlay_video
                    writeVideo(overlayV, overlayImage);
                end
                if config.video.overlay_frames
                    imwrite(overlayImage, fullfile(overlayFolder, expected(i)));
                end
            end
        end
        % Close output videos.
        if ~isempty(V)
            close(V);
            V = [];
        end
        if ~isempty(overlayV)
            close(overlayV);
            overlayV = [];
        end

        %% CONVERT RULER COORDINATES TO DEGREES
        % Retain the original ruler projection, angle conversion and zero convention.
        % The laser reflection doubles the mirror deflection, hence the 0.5.
        % Subtract the reference angle before converting radians to degrees.
        tan1 = atan2((position - xC), D);
        tan2 = atan2((x0 - xC), D);
        theta_all = 0.5*(tan1 - tan2)*180/pi; % convert to degrees
        % Collect every affected 1-based frame number rather than stopping at
        % the first missing laser. No interpolation or partial CSV is accepted.
        metadata.invalid_frames = find(~isfinite(position) | ~isfinite(theta_all))';
        assert(isempty(metadata.invalid_frames), ...
            'Missing laser/nonfinite angles at frames %s', mat2str(metadata.invalid_frames));
        %% AUDIO TIMING AND INDIVIDUAL DATA EXPORT
        mov_name = fullfile(config.paths.trials, file + config.video.extension);
        % Save fresh audio candidates separately from the preserved manual
        % corrections. Only the effective startup time shifts comparison plots.
        metadata.audio_candidates_s = estimateSoundTimes(mov_name, config.audio, [t, theta_all]);
        startup = effectiveTime(metadata, 'startup');
        effectiveTime(metadata, 'shutdown'); % Validate both events before accepting the trial.
        % Write theta data to a CSV file, keeping the original column format.
        dataFinal = [t, theta_all];
        outputTable = array2table(dataFinal, ...
            'VariableNames', {'time_seconds', 'angle_degrees'});
        writetable(outputTable, csvPath);
        % Publish validity after the CSV write succeeds; failures below remove
        % that CSV again and leave metadata invalid for Python to reject.
        metadata.valid = true;
        metadata.failure_reason = '';
        writeMetadata(metadataPath, metadata);
        validTrials(end+1) = struct('id', file, 'group', groups(trialIndex), ...
            'time', t, 'angle', theta_all, 'startup', startup); %#ok<SAGROW>
        fprintf('Saved %s\n', csvPath);
    catch err
        if ~isempty(V)
            close(V);
        end
        if ~isempty(overlayV)
            close(overlayV);
        end
        if isfile(csvPath)
            delete(csvPath);
        end
        % A malformed existing JSON is not overwritten: keep manual edits recoverable.
        if isfield(metadata, 'valid')
            metadata.valid = false;
            metadata.failure_reason = err.message;
            writeMetadata(metadataPath, metadata);
        end
        fprintf(2, 'FAILED %s: %s\n', file, err.message);
    end
end

%% PLOT INDIVIDUAL TRIALS
% One line per trial, on its own original samples; audio shift is display-only.
for group = unique(groups(:))'
    fig = figure;
    hold on;
    lines = 0;
    for i = 1:numel(validTrials)
        trial = validTrials(i);
        if trial.group == group
            plot(trial.time - trial.startup, trial.angle, ...
                'LineWidth', 1.5, 'DisplayName', trial.id);
            lines = lines + 1;
        end
    end
    if lines == 0
        close(fig);
        continue;
    end
    xlabel('Time from effective startup (s)');
    ylabel('Angular displacement (degrees)');
    title(group + ": individual angular displacements");
    legend('show');
    grid on;
    exportgraphics(fig, fullfile(config.paths.output, group, 'angular_displacement.png'));
end

function writeMetadata(path, metadata)
% Write the one shared per-trial JSON, including any preserved user selections.
% onCleanup closes the file even if encoding or writing raises an error.
    fid = fopen(path, 'w');
    assert(fid ~= -1, 'Cannot write metadata: %s', path);
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', jsonencode(metadata, 'PrettyPrint', true));
end

function value = effectiveTime(metadata, event)
% MATLAB counterpart of Python's effective_event_time: prefer a manual cue,
% otherwise use the candidate. Both are original video seconds, independent
% of the torque integration bounds stored elsewhere in this same metadata.
    value = metadata.audio_candidates_s.(event);
    if isfield(metadata.manual_event_times_s, event) && ~isempty(metadata.manual_event_times_s.(event))
        value = metadata.manual_event_times_s.(event);
    end
    validateattributes(value, {'numeric'}, {'scalar', 'real', 'finite'});
end

function outputImage = drawHollowBlackCircle(inputImage, centerX, centerY, radius, thickness)
% Diagnostic overlay only; the original tracking image is not modified.
    outputImage = inputImage;
    [imageHeight, imageWidth, ~] = size(outputImage);

    [gridX, gridY] = meshgrid(1:imageWidth, 1:imageHeight);
    distanceFromCenter = sqrt((gridX - centerX).^2 + (gridY - centerY).^2);
    ringMask = distanceFromCenter >= (radius - thickness) & distanceFromCenter <= radius;

    if isinteger(outputImage)
        blackValue = intmax(class(outputImage));
        blackValue = blackValue - blackValue;
    else
        blackValue = 0;
    end

    for channel = 1:size(outputImage, 3)
        channelImage = outputImage(:,:,channel);
        channelImage(ringMask) = blackValue;
        outputImage(:,:,channel) = channelImage;
    end
end
