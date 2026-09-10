clc; clear; clf; close all;

trialCount = 1;
file = "1000r";

motionThreshold = 2;   % degrees / second

%% READ ALL TRIALS

trialData = cell(1, trialCount);

for trial = 1:trialCount
    name = "data/" + file + int2str(trial) + "_data.csv";
    trialData{trial} = readmatrix(name);
end

% Get dt from the actual saved timestamps
dt = median(diff(trialData{1}(:,1)));

fprintf("dt = %.6f s\n", dt);


%% FIND WHEN MOTION STARTS IN EACH TRIAL

devLoc = zeros(1, trialCount);
devTime = zeros(1, trialCount);

for trial = 1:trialCount

    tTrial = trialData{trial}(:,1);
    thetaTrial = trialData{trial}(:,2);

    % Angular velocity estimate
    deriv = gradient(thetaTrial, tTrial);

    % Use ABSOLUTE derivative so this works for forward and reverse
    idx = find(abs(deriv) > motionThreshold, 1, 'first');

    if isempty(idx)
        error("No motion detected in trial %d", trial);
    end

    devLoc(trial) = idx;
    devTime(trial) = tTrial(idx);
end

disp("Detected motion start times:")
disp(devTime)


%% ALIGN TRIALS

% Earliest detected start
earliestStart = min(devTime);

% Amount each trial should be shifted left
devDiff = devTime - earliestStart;

fprintf("Time shifts:\n")
disp(devDiff)


%% DETERMINE COMMON TIME DOMAIN

% Find the latest time available after shifting for each trial
endTimes = zeros(1, trialCount);

for trial = 1:trialCount
    tTrial = trialData{trial}(:,1);
    tShift = tTrial - devDiff(trial);

    endTimes(trial) = max(tShift);
end

% Only use times for which EVERY trial has data
commonEnd = min(endTimes);

tCommon = (0:dt:commonEnd)';


%% INTERPOLATE EACH TRIAL ONTO COMMON TIME GRID

interpData = NaN(length(tCommon), trialCount);

for trial = 1:trialCount

    tTrial = trialData{trial}(:,1);
    thetaTrial = trialData{trial}(:,2);

    tShift = tTrial - devDiff(trial);

    % Ignore any NaNs that may occur in the raw data
    valid = isfinite(tShift) & isfinite(thetaTrial);

    interpData(:,trial) = interp1( ...
        tShift(valid), ...
        thetaTrial(valid), ...
        tCommon, ...
        'pchip', ...
        NaN);
end


%% KEEP ONLY POINTS VALID FOR ALL THREE TRIALS

allValid = all(isfinite(interpData), 2);

tCommon = tCommon(allValid);
interpData = interpData(allValid,:);


%% COMPUTE MEAN

meanData = mean(interpData, 2);


%% PLOT

f = figure(1);
f.Theme = "Light";

hold on

% Plot shifted raw trials
for trial = 1:trialCount

    tTrial = trialData{trial}(:,1);
    thetaTrial = trialData{trial}(:,2);

    tShift = tTrial - devDiff(trial);

    plot( ...
        tShift, ...
        thetaTrial, ...
        '-.', ...
        'LineWidth', 1.5, ...
        'DisplayName', "Trial " + int2str(trial));
end

% Plot mean
plot( ...
    tCommon, ...
    meanData, ...
    '-', ...
    'LineWidth', 2.5, ...
    'DisplayName', "Mean of All Trials");

hold off

legend()
xlabel('Time (s)');
ylabel('Angular Displacement (degrees)');
title('Trial Data Comparison');


%% SAVE MEAN DATA

name = "data/" + file + "_meanData.csv";

dataFinal = [tCommon, meanData];

writematrix(dataFinal, name, "Delimiter", ',');
type(name);