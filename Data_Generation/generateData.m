
clc
clf
close all 
clear

for trial = [1, 2, 3]
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
    % 
    % if ~isfolder(saveFolder)
    %     mkdir(saveFolder);
    % end
    % 
    % i = 1;
    % 
    % while hasFrame(v)
    %     img = readFrame(v);
    %     filename = sprintf("%03d",i)+".jpg";
    %     fullname = file + "_frames/" + filename;
    %     imwrite(img,fullname)    % Write to a JPEG file (001.jpg, 002.jpg, ..., 121.jpg)
    %     i = i+1;
    % end
  
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

    % Video Initialization
    vidName = file + "_maskedData.mp4";
    V = VideoWriter(vidName, 'MPEG-4');
    open(V)

    while isfile(path)
        % Reads jpg as image matrix
        img = imread(path);
        % Masks image using mask defined in Color Thresholder app. This app
        % is in the image processing toolbox for matlab.
        % Outputs a logical array (BW) and the masked image in color (matrix) 
        [BW,maskedImage] = ycbcrMask(img); 

        % Writes masked image to a video 
        writeVideo(V, maskedImage)

        % This flattens the masked image from 3 dimensional to 2
        % dimensional. We don't care where the color is on the spectrum
        % after the mask, we simply care about its magnitude.
        flattenImage = sum(maskedImage, 3); 

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

        % Flattens the logical array BW to the x dimension
        xBW = mean(BW, 1);

        % Same as above, just for BW image
        [row, col] = find(xBW);

        % Takes means of BW image wrt axis, but for x position, we use the
        % magnitude of color at the position as a weight for a weighted
        % mean - the laser is more likely to be somewhere there is higher
        % magnitude of color than not. 
        yMid = mean(row);
        xMid = mean(col, 'Weights', weight);

        % If first iteration, define zero as the mean x position of laser
        if findZero
            zero = xMid;
            findZero = false;
        end

        % Apped current position to data array, normalized by zero position
        position = [position; (xMid - zero)];

        % Define current time as i*dt
        t = [t; (i-1)*dt];

        % Update i, filename, and path before next loop
        i = i + 1;
        filename = sprintf("%03d",i)+".jpg";
        path = cd + "/" + folder + filename;
    end

    % Close maskedImage video
    close(V)


    %%

    % Experimental setup hard-coded values from Kelly
    fixed_distance = 49; % Distance from sprinkler mirror to ruler board
    inches_per_pixel = 28/1920; % Number of inches in frame over horizontal resolution


    % Convert raw data from iamge to inches - position is the pixel index
    inches = position*inches_per_pixel;% Convert horizontal pixel position to physical displacement (inches)

    % Convert inches to degrees
    theta_all = 0.5*atan2(inches, fixed_distance)*180/pi; % For whole run


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