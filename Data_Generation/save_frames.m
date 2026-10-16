
clc
clf
close all 
clear

% Resolve files relative to this script, independently of MATLAB's current folder.
scriptDir = fileparts(mfilename('fullpath'));
trialsDir = fullfile(scriptDir, "trials");
framesDir = fullfile(scriptDir, "frames");
calibrationDir = fullfile(scriptDir, "calibration");
calibrationFrame = fullfile(calibrationDir, "frame.jpg");

% Extract each batch's frames and refresh its calibration reference image.
% Frame extraction takes time and only needs to be done once per batch.

% The first listed trial supplies the calibration image for this batch.
trials = [1];
for trialIndex = 1:numel(trials)
    trial = trials(trialIndex);
    file = "1000f" + int2str(trial);
    
    %% GET FRAMES
    
    mov_name = fullfile(trialsDir, file + ".MTS");
    
    if ~isfile(mov_name)
        error("Recording not found: " + mov_name)
    end
    
    % create video reader object for reading the video files
    v = VideoReader(mov_name); 
    if ~hasFrame(v)
        error("Recording contains no readable frames: " + mov_name)
    end
 
    saveFolder = fullfile(framesDir, file + "_frames");
    
    if ~isfolder(saveFolder)
        mkdir(saveFolder);
    end
    
    i = 1;
    
    while hasFrame(v)
        img = readFrame(v);
        if trialIndex == 1 && i == 1
            if ~isfolder(calibrationDir)
                mkdir(calibrationDir);
            end
            imwrite(img, calibrationFrame);
            fprintf("Calibration frame saved from %s to %s\n", mov_name, calibrationFrame);
        end
        filename = sprintf("%03d",i)+".jpg";
        fullname = fullfile(saveFolder, filename);
        imwrite(img,fullname)    % Write to a JPEG file (001.jpg, 002.jpg, ..., 121.jpg)
        i = i+1;
    end
end
