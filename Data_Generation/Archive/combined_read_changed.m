%plz specify a t=0 (when the flow turns on) and a t>t0 where the spring motion has mostly decayed. These times are verified by eye using the video. \
% define a function to get the video without returning anything
function [] = combined_read_changed(input_filename)
% read the video
v = VideoReader(input_filename);
fps = v.FrameRate;
loc = zeros(v.NumFrames, 2); % Empty array to store tracked locations 
dt = 1/fps; % Time step size
num_frames = v.NumFrames;
top_cutoff = 0.7; % Fraction of the frame from the top for code to look at
bottom_cutoff = 0.9; % Fraction of the frame from the bottom for code to look at


% Display first frame and ask the user to draw a box covering the lazer dot
frame = readFrame(v); % Read first frame
r = frame(1:(bottom_cutoff*v.Height),:, 1);%:, :, 1); % Get red channel only with the frame we cut
imshow(r); % Show first frame
temp = getrect; % Get rectangle
marker_diameter = floor(temp(4)); 
% Use the height of the bounding box as the diameter of the laser dot (rounded down)
%can change to use width by changing to temp(3)
v.CurrentTime = 0; % Reset to read from time = 0



% Perform tracking
i = 1;
tStart = tic; 
while hasFrame(v)
    %if there is still frame not being read then we read the next frame
    frame = readFrame(v); 
    % Get red channel only with the frame we cut
    r = frame(1:(bottom_cutoff*v.Height), :, 1);
    %r = frame((top_cutoff*v.Height):(bottom_cutoff*v.Height), :, 1);
    detect = conv2(r, ones(marker_diameter), 'same'); % 2D moving sum
    %[max_value, p] = max(detect(:)); % maxvalue is the value and p is index

    % this is the reason why this may not be the middle of the blur
    % the max function will retrun the first max value if there are more than one max

    %change our index to the 2d coordinate since detect(:) make a column vector
    %[py, px] = ind2sub(size(detect), p); 
    %loc(i,:) = [py, px]; % Store peak location




    % but if we want to get the middle point when it's a blur I changed the code as
    threshold = 0.95 * max(detect(:));
    [rows, cols] = find(detect > threshold);
    weights = detect(detect > threshold);
    py = sum(rows .* weights) / sum(weights);
    px = sum(cols .* weights) / sum(weights);
    loc(i,:) = [py, px];



    
    % This num↓ indicates every 'n' frames show progress in command window. Make very large if you want to turn off progress updates.
    if mod(i, 1) == 0
         %If you want to check the current progress
        sprintf("Tracking... current progress: %d/%d", i, v.NumFrames)
        imshow(r); 
        hold on
        rectangle('Position', [px-floor(marker_diameter/2), ...
                               py-floor(marker_diameter/2), ...
                               marker_diameter, marker_diameter], ... 
                  'EdgeColor', 'g', 'LineWidth', 2); 
        drawnow; 
        pause(1);

    end
    i = i+1; 
end
tEnd = toc(tStart); 
clf


    start = round(input("Please input time at which flow started/stopped: ")*fps)
    stop = round(input("Please input a time at which transient state has mostly decayed: ")*fps)
    fixed_distance = 49; % Distance from sprinkler mirror to ruler board
    inches_per_pixel = 28/1920; % Number of inches in frame over horizontal resolution
    





    % Defining deflection angle and time for plotting
    inches = loc(:,2)*inches_per_pixel;% Convert horizontal pixel position to physical displacement (inches)
    theta_all = 0.5*atan2(inches, fixed_distance); % For whole run
    theta_raw = theta_all(start:stop); % For selected time window
    %theta_asymptote = mean(theta_all(round(stop-10*fps):stop)); % Neutral point of oscillator

    % if just use a short video
    start_idx = max(1, round(stop - 10 * fps)); 
    theta_asymptote = mean(theta_all(start_idx:stop)); % Neutral point of oscillator

    theta = theta_raw-theta_asymptote; % Removes offset from starting location and zeroes to asymptote
    time = transpose(linspace(0, (stop-start+1)/fps, stop-start+1)); % Start time at 

    
    
    
    
    % Get local max/min (Separation is important to prevent noise from being considered a max/min. You may need to tweak sensitivity.)
    TF_max = islocalmax(theta,'MinSeparation', 0.5 ,'SamplePoints',time); % Get all peaks
    TF_min = islocalmin(theta,'MinSeparation', 0.5 ,'SamplePoints',time); % Get all valleys
    TF = or(TF_min, TF_max); % Saves incides of peaks and valleys

    % Plot deflection angle of the sprinkler over desired time and circle all peaks and valleys
    hold on;
    grid on;
    plot(time, theta*180/pi, time(TF), theta(TF)*180/pi,'o', Linewidth = 1.5)
    title("Spring Response")
    xlabel('Time Since Pump Start (s)')
    ylabel('Angular Position (deg)')

    % Use absolute value of peaks and valleys to fit a bounding exponential decay
    fo = fitoptions('Method', 'NonlinearLeastSquares');
    myfittype = fittype("a*exp(-b*time)",...
    dependent="theta",independent="time",...
    coefficients=["a" "b"])
    curve = fit(time(TF), abs(theta(TF))*180/pi, myfittype) % Can include start points here for the coefficients
    
    % Plot all peaks and valleys (absolute value) and fitted exponential
    figure
    hold on;
    grid on;
    plot(curve, time(TF), abs(theta(TF))*180/pi)
    title("Fitted Decay: Standard Sprinkler, Re = XXX")
    xlabel('Time Since Pump Start (s)')
    ylabel('Angular Position (deg)')
    hold off;
    theta_radians = theta;
    

    T = table(time, theta_radians, 'VariableNames', {'time (sec)', 'theta (rad)'});
    writetable(T, "rev_1500_trial2.csv");

end
