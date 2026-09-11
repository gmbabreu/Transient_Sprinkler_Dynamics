
clc
clf
close all 
clear

D = 48.5; % Fixed Distance from Camera to mirror at shortest point
saveData = true;

overlayVideo = true;  % Save the overlayed video
    overlayFrame = rand < 0;  % Save individual overlayed frames, default off
circleRadius = 8;  
circleThickness = 2;

%% Calibration
calibrationFile = "frame_rulerCalibration.mat";

trialCount = 1;

for trial = 1:1:trialCount
    file = "1000f" + int2str(trial);
    
    % Get movie data
    mov_name = "trials/" + file + ".MTS";
    
    if ~isfile(mov_name)
        error(mov_name + " is not in working directory. Please move move or code to working directory!")
    end

    % create video reader object for reading the video files
    v = VideoReader(mov_name); 
    fps = v.FrameRate;
    dt = 1/fps; % Time step size
    fprintf("Time step size is %.4f seconds\n", dt);

    if ~isfile(calibrationFile)
        error("Missing calibration file " + calibrationFile + ". Run calibrateRulerAxis.m first.")
    end

    calibrationData = load(calibrationFile, "rulerCalibration");
    rulerCalibration = calibrationData.rulerCalibration;
    xC = rulerCalibration.originValue;   % 18 inches
  
    %% READ, MASK, NAD EXTRACT DATA

    folder = fullfile("frames", file + "_frames");

    pink = [170/225, 22/225, 181/225]; % Some colors I like
    darkpink = [137/225, 18/225, 145/225];

    % Defining loops for names to read through frames
    i = 1;
    filename = sprintf("%03d",i)+".jpg";
    path = fullfile(pwd, folder, filename);

    if ~isfile(path)
        error("First frame not found: " + path + ". Run save_frames.m for this trial first.");
    end

    % This will enter a loop to define the starting x-position as "0"
    findZero = true;

    % Initializing arrays for position of laser (pos) and time (t)
    position = [];
    t = [];

    % Always save the masked video as before.
    if ~isfolder("data")
        mkdir("data");
    end
    vidName = "data/" + file + "_maskedData.mp4";
    V = VideoWriter(vidName, 'MPEG-4');
    V.FrameRate = fps;
    V.Quality = 95;
    open(V)

    if overlayVideo
        overlayVidName = "data/" + file + "_overlayData.mp4";
        overlayV = VideoWriter(overlayVidName, 'MPEG-4');
        overlayV.FrameRate = fps;
        overlayV.Quality = 95;
        open(overlayV)
    end

    if overlayFrame
        overlayFolder = "data/" + file + "_overlayFrames/";
        if ~isfolder(overlayFolder)
            mkdir(overlayFolder);
        end
    end

    while isfile(path)
        % Reads jpg as image matrix
        img = imread(path);
        % Masks image using mask defined in Color Thresholder app. This app
        % is in the image processing toolbox for matlab.
        % Outputs a logical array (BW) and the masked image in color (matrix) 
        [BW,maskedImage] = ycbcrMask(img); 

        % Find the predicted laser position from the masked pixels.
        laserImage = sum(maskedImage, 3);
        [laserRows, laserCols] = find(BW);
        laserWeights = double(laserImage(BW));

        if isempty(laserCols)
            xMid = NaN;
            yMid = NaN;
        else
            xMid = mean(laserCols, 'Weights', laserWeights);
            yMid = mean(laserRows, 'Weights', laserWeights);
        end

        rulerInches = projectToRulerAxis([xMid, yMid], rulerCalibration);

        % Overlay a hollow black circle on the original frame at the
        % predicted laser location.
        overlayImage = img;
        if ~isnan(xMid) && ~isnan(yMid)
            overlayImage = drawHollowBlackCircle(overlayImage, xMid, yMid, circleRadius, circleThickness);
        end

        writeVideo(V, maskedImage)

        if overlayVideo
            writeVideo(overlayV, overlayImage)
        end

        if overlayFrame
            frameName = sprintf("%03d", i) + ".jpg";
            imwrite(overlayImage, overlayFolder + frameName);
        end

        % This flattens the masked image from 3 dimensional to 2
        % dimensional. We don't care where the color is on the spectrum
        % after the mask, we simply care about its magnitude.
        flattenImage = laserImage; 

        % Now we flatten the image in the y dimension, as we only care
        % about its position along the x axis
        xFlat = mean(flattenImage, 1);

        % Outputs all the row indicies (xRow) and column indicies where
        % xFlat is nonzero ( remember it is masked - so only nonzero where
        % laser is )
        [xRow, xCol] = find(xFlat);

        % Defines a new vector weight, which consists of xFlat at every
        % index defined in xCol.
        weight = xFlat(xCol);

        % Takes means of BW image wrt axis, but for x position, we use the
        % magnitude of color at the position as a weight for a weighted
        % mean - the laser is more likely to be somewhere there is higher
        % magnitude of color than not. 
        % xMid and yMid are already computed from the full mask above.

        % If first iteration, define zero as the ruler position of the laser
        if findZero && isfinite(rulerInches)
            x0 = rulerInches;
            findZero = false;
        end

        if findZero
            error("Laser was never detected; could not define x0.")
        end

        % Apped current position to data array
        position = [position; rulerInches];

        % Define current time as i*dt
        t = [t; (i-1)*dt];

        % Update i, filename, and path before next loop
        i = i + 1;
        filename = sprintf("%03d",i)+".jpg";
        path = fullfile(folder, filename);
    end

    % Close output videos
    close(V)
    if overlayVideo
        close(overlayV)
    end


    %%

    % % OLD OLD OLD OLD
    % % Experimental setup hard-coded values from Kelly
    % fixed_distance = 49; % Distance from sprinkler mirror to ruler board
    % inches_per_pixel = 28/1920; % Number of inches in frame over horizontal resolution


    % Convert ruler coordinates in inches to degrees.
    tan1 = atan2((position - xC), D);
    tan2 = atan2((x0 - xC), D);
    theta_all = 0.5*(tan1 - tan2)*180/pi; % convert to degrees


    if saveData
        % Set up to write theta data to a txt file 
        name = fullfile("data", file + "_data.csv");
        dataFinal = [t, theta_all];
        blocker = false;

        % This loop ensures user does not overwrite data. Does require user to
        % manual say yes/no
        if isfile(name)
            disp("Previous data will be deleted, would you like to proceed? [yes/no]")
            prompt = "=> ";
            response = input(prompt, 's');
            if response == "yes"
                blocker = true;
            elseif response == "no"
                    error("Please delete or move old data. Operation aborted by user")
            else
                prompt = "Please provide [yes/no], this is case and spacing senesitive";
                response = input(prompt, 's');
                if response == "yes"
                    blocker = true;
                elseif response == "no"
                        error("Please delete or move old data. User has ended the code.")
                end 
            end 
        end 

        if blocker
            delete(name); % Delete old data if user confirmed
        end

        % write theta data to a txt file 
        outputTable = array2table(dataFinal, ...
            'VariableNames', {'time_seconds', 'angle_degrees'});
        writetable(outputTable, name);
        type(name);
    end
    

end

function outputImage = drawHollowBlackCircle(inputImage, centerX, centerY, radius, thickness)
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