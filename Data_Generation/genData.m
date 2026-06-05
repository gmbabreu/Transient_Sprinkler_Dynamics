
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
    
    % create video reader object for reading teh video files
    v = VideoReader(mov_name); 
    fps = v.FrameRate;
    dt = 1/fps; % Time step size
    
    %% READ, MASK, NAD EXTRACT DATA
    
    % Display first frame and ask the user to draw a box covering the lazer dot
    frame = readFrame(v); % Read first frame
    
    saveFolder = file + "_frames/";
    
    if ~isfolder(saveFolder)
        mkdir(saveFolder);
    end
    
    i = 1;
    
    while hasFrame(v)
        img = readFrame(v);
        filename = sprintf("%03d",i)+".jpg";
        fullname = file + "_frames/" + filename;
        imwrite(img,fullname)    % Write to a JPEG file (001.jpg, 002.jpg, ..., 121.jpg)
        i = i+1;
    end
    
    %% READ FRAMES
    
    folder = file + "_frames/";
    
    pink = [170/225, 22/225, 181/225];
    darkpink = [137/225, 18/225, 145/225];
    
    i = 1;
    filename = sprintf("%03d",i)+".jpg";
    path = cd + "/" + folder + filename;
    
    % savefolder = cd + "/" + file + "_masked/";
    
    
    
    
    findZero = true;
    
    
    pos = [];
    t = [];
    
    vidName = file + "_maskedData.mp4";
    
    V = VideoWriter(vidName, 'MPEG-4');
    open(V)
    
    while isfile(path)
        img = imread(path);
        [BW,maskedRGBImage] = laserMask(img); % masks photo by the red color (user defined)
    
    
        [row, col] = find(BW);
        % imwrite(maskedRGBImage, saveName);
        writeVideo(V, maskedRGBImage)
    
    
        yMid = mean(row);
        xMid = mean(col);
        
        if findZero
            zero = xMid;
            findZero = false;
        end
    
    
        pos = [pos; (xMid - zero)];
        t = [t; i*dt];
    
    
        i = i + 1;
        filename = sprintf("%03d",i)+".jpg";
        path = cd + "/" + folder + filename;
    end
    
    close(V)
    
    
    %%
    
    
    fixed_distance = 49; % Distance from sprinkler mirror to ruler board
    inches_per_pixel = 28/1920; % Number of inches in frame over horizontal resolution
    
    
    
    % Defining deflection angle and time for plotting
    inches = pos*inches_per_pixel;% Convert horizontal pixel position to physical displacement (inches)
    theta_all = 0.5*atan2(inches, fixed_distance)*180/pi; % For whole run
    
    
    
    name = file + "_data.txt";
    
    dataFinal = [t, theta_all];
    
    blocker = false;
    
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
    
    
    writematrix(dataFinal, name, "Delimiter", '\t');
    type(name);
    
    % 
    % fig = figure(1);
    % fig.Theme = "Light";
    % 
    % plot(t, theta_all, "-o", "Color", pink, 'LineWidth', 3, 'MarkerEdgeColor', darkpink, 'MarkerSize', 3)
end