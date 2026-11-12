% Extract the configured groups/repeats before calibration/data generation.
% loadBatchConfig expands trial_groups and trials_per_group into full trial IDs.
% Each folder contains numbered JPEGs plus a .complete marker with frame count
% and frame rate. Set video.regenerate_frames in experiment_config.json when
% a recording changes; the marker is a completion check, not a content hash.
clc
close all
clear

% Resolve files relative to this script, independently of MATLAB's current folder.
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
config = loadBatchConfig();

for trialIndex = 1:numel(config.trials)
    file = config.trials(trialIndex);

    %% GET FRAMES
    % Frame extraction takes time and only needs to be done once per trial.
    mov_name = fullfile(config.paths.trials, file + config.video.extension);
    saveFolder = fullfile(config.paths.frames, file + "_frames");
    marker = fullfile(saveFolder, '.complete');
    try
        complete = false;
        % A marker alone is insufficient: require exactly the expected JPEGs
        % as well. Missing or extra frames cause the folder to be rebuilt.
        if isfile(marker)
            info = sscanf(fileread(marker), '%f'); % frame count, frame rate
            if numel(info) == 2 && all(isfinite(info)) && info(1) > 0 && info(1) == fix(info(1)) && info(2) > 0
                expected = compose('%03d.jpg', (1:info(1))');
                complete = numel(dir(fullfile(saveFolder, '*.jpg'))) == info(1) ...
                    && all(isfile(fullfile(saveFolder, expected)));
            end
        end
        if complete && ~config.video.regenerate_frames
            fprintf('Skipping complete extraction: %s\n', file);
        else
            assert(isfile(mov_name), 'Recording not found: %s', mov_name);
            % A changed extraction requires fresh tracking before Python can use it.
            % Preserve manual timing/selection fields while invalidating the
            % metadata and removing angles computed from the previous frames.
            metadataPath = fullfile(config.paths.data, file + "_metadata.json");
            if isfile(metadataPath)
                metadata = jsondecode(fileread(metadataPath));
                metadata.valid = false;
                metadata.failure_reason = 'Frames changed; rerun generateData.m';
                fid = fopen(metadataPath, 'w');
                assert(fid ~= -1, 'Cannot invalidate metadata: %s', metadataPath);
                fprintf(fid, '%s\n', jsonencode(metadata, 'PrettyPrint', true));
                fclose(fid);
            end
            csvPath = fullfile(config.paths.data, file + "_data.csv");
            if isfile(csvPath)
                delete(csvPath);
            end
            if ~isfolder(saveFolder)
                mkdir(saveFolder);
            end
            % Clear only generated frames and the marker before ANY rewrite.
            if isfile(marker)
                delete(marker);
            end
            oldFrames = dir(fullfile(saveFolder, '*.jpg'));
            for j = 1:numel(oldFrames)
                delete(fullfile(saveFolder, oldFrames(j).name));
            end
            % Create a video reader object for reading the video file.
            v = VideoReader(mov_name);
            count = 0;
            % Read sequentially to exhaustion instead of relying on an estimated
            % frame count from duration. Names retain the original 001.jpg style.
            while hasFrame(v)
                img = readFrame(v);
                count = count + 1;
                imwrite(img, fullfile(saveFolder, sprintf('%03d.jpg', count)));
            end
            assert(count > 0, 'Recording has no readable frames: %s', mov_name);
            % Write the marker last. Interrupted extraction leaves no marker,
            % so generateData will refuse to treat a partial folder as complete.
            fid = fopen(marker, 'w');
            assert(fid ~= -1, 'Cannot write completion marker: %s', marker);
            fprintf(fid, '%d %.17g\n', count, v.FrameRate);
            fclose(fid);
            fprintf('Extracted %d frames: %s\n', count, file);
        end
        % The first trial supplies the batch reference, including on a skipped run.
        % Never silently replace a reference used by an existing ruler calibration.
        if trialIndex == 1 && ~isfile(config.paths.calibration_reference)
            referenceDir = fileparts(config.paths.calibration_reference);
            if ~isfolder(referenceDir)
                mkdir(referenceDir);
            end
            copyfile(fullfile(saveFolder, '001.jpg'), config.paths.calibration_reference);
        end
    catch err
        % A missing/unreadable recording must not prevent later trials running.
        % Do not substitute a later trial for the first trial's batch reference.
        fprintf(2, 'FAILED %s: %s\n', file, err.message);
    end
end
