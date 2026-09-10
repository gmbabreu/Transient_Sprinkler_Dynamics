
clc
clf
close all 
clear

% THE FOLLOWING MAKES FRAME DIRECTORIES AND FRAMES
% This only needs to be done once on your personal device. It is necessary
% for this code to run to have frames. This process takes time.

for trial = [1]
    file = "1000f" + int2str(trial);
    
    %% GET FRAMES
    
    mov_name = "trials/" + file + ".MTS";
    
    if ~isfile(mov_name)
        error(mov_name + " is not in working directory. Please move move or code to working directory!")
    end
    
    % create video reader object for reading the video files
    v = VideoReader(mov_name); 
 
    saveFolder = fullfile("frames", file + "_frames");
    
    if ~isfolder(saveFolder)
        mkdir(saveFolder);
    end
    
    i = 1;
    
    while hasFrame(v)
        img = readFrame(v);
        filename = sprintf("%03d",i)+".jpg";
        fullname = fullfile(saveFolder, filename);
        imwrite(img,fullname)    % Write to a JPEG file (001.jpg, 002.jpg, ..., 121.jpg)
        i = i+1;
    end
end