
clc
clf
close all 
clear

D = 50; % Fixed Distance from Camera to mirror at shortest point
xC = 40000; % x position on ruler at shortest point (point where D is measured)
    % xC is in inches and must match the calibration-defined ruler units.
saveData = true;

overlayVideo = true;  % Save the overlayed video
    overlayFrame = rand < 0;  % Save individual overlayed frames, default off
circleRadius = 8;  
circleThickness = 2;


for trial = [1, ]
    file = "1500f" + int2str(trial);
    
    %% GET FRAMES
    
    mov_name = file + ".mov";
    
    if ~isfile(mov_name)
        error(mov_name + " is not in working directory. Please move move or code to working directory!")
    end
    
    % create video reader object for reading the video files
    v = VideoReader(mov_name); 
    fps = v.FrameRate;
    dt = 1/fps; % Time step size
    
%%%%%%%%%%%%%%%%%%%%%%% I M P O R T A N T %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE FOLLOWING MAKES FRAME DIRECTORIES AND FRAMES
% This only needs to be done once on your personal device. It is necessary
% for this code to run to have frames. By uncommenting below, you will make
% frames for trail as it runs. This process takes time.

% Simply uncomment lines 31 through 45 to make directories and run code as
% normal. 

 
    % saveFolder = file + "_frames/";
    
    % if ~isfolder(saveFolder)
    %     mkdir(saveFolder);
    % end
    
    % i = 1;
    
    % while hasFrame(v)
    %     img = readFrame(v);
    %     filename = sprintf("%03d",i)+".jpg";
    %     fullname = file + "_frames/" + filename;
    %     imwrite(img,fullname)    % Write to a JPEG file (001.jpg, 002.jpg, ..., 121.jpg)
    %     i = i+1;
    % end

    %% Calibration
    calibrationFile = file + "_rulerCalibration.mat";

    if ~isfile(calibrationFile)
        error("Missing calibration file " + calibrationFile + ". Run calibrateRulerAxis.m first.")
    end

    calibrationData = load(calibrationFile, "rulerCalibration");
    rulerCalibration = calibrationData.rulerCalibration;
  
    %% READ, MASK, NAD EXTRACT DATA

    folder = file + "_frames/";

    pink = [170/225, 22/225, 181/225]; % Some colors I like
    darkpink = [137/225, 18/225, 145/225];

    % Defining loops for names to read through frames
    i = 1;
    filename = sprintf("%03d",i)+".jpg";
    path = cd + "/" + folder + filename;

    % This will enter a loop to define the starting x-position as "0"
    findZero = true;

    % Initializing arrays for position of laser (pos) and time (t)
    position = [];
    t = [];

    % Always save the masked video as before.
    vidName = file + "_maskedData.avi";
    V = VideoWriter(vidName, 'Motion JPEG AVI');
    open(V)

    if overlayVideo
        overlayVidName = file + "_overlayData.avi";
        overlayV = VideoWriter(overlayVidName, 'Motion JPEG AVI');
        open(overlayV)
    end

    if overlayFrame
        overlayFolder = file + "_overlayFrames/";
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
            xMid = round(mean(laserCols, 'Weights', laserWeights));
            yMid = round(mean(laserRows, 'Weights', laserWeights));
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
        if findZero
            x0 = rulerInches;
            findZero = false;
        end

        % Apped current position to data array
        position = [position; rulerInches];

        % Define current time as i*dt
        t = [t; (i-1)*dt];

        % Update i, filename, and path before next loop
        i = i + 1;
        filename = sprintf("%03d",i)+".jpg";
        path = cd + "/" + folder + filename;
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
    theta_all = 0.5*(tan1 - tan2)*180/pi; % For whole run


    if saveData
        % Set up to write theta data to a txt file 
        name = file + "_data.txt";
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
        writematrix(dataFinal, name, "Delimiter", '\t');
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