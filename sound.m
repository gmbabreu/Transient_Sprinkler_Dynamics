% Inspect audio cues for pump startup/shutdown, in original recording seconds.
% A loudest sample is only a candidate cue, NOT a motor-onset detector.
% First run: view the overview. Then narrow startTime/endTime around one
% event and enable playback. Repeat for the other event.
repoDir = fileparts(mfilename('fullpath'));
trialName = '1000r1';
videoFileName = fullfile(repoDir, 'Data_Generation', 'trials', [trialName '.MTS']);
startTime = 9.5;
endTime = 13;             % Inf means the end of the audio
playbackEnabled = false;   % enable after narrowing the inspection window
playbackChannels = [1 2];  % audition other channels separately if needed
playbackGain = 0.25;       % fixed attenuation; start with low speaker volume
audioToVideoOffset_s = 0;  % video time = audio time + offset; verify sync

% Read every channel for analysis, independently of playback selection.
disp(['Reading file: ', videoFileName, ' ...']);
[audioData, fs] = audioread(videoFileName);
numChannels = size(audioData, 2);
disp(['Number of audio channels: ', num2str(numChannels)]);
totalSamples = size(audioData, 1);
duration_s = totalSamples / fs;
endTime = min(endTime, duration_s);
validateattributes(startTime, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
validateattributes(endTime, {'numeric'}, {'scalar', 'real', 'finite', '>', startTime});
assert(startTime < duration_s, 'startTime must be inside the recording.');

% Half-open interval [startTime, endTime), using actual sample timestamps.
startIndex = ceil(startTime * fs) + 1;
endIndex = min(totalSamples, ceil(endTime * fs));
assert(endIndex >= startIndex, 'Inspection window contains no audio samples.');
audioSegment = audioData(startIndex:endIndex, :);
timeVector = ((startIndex:endIndex)' - 1) / fs + audioToVideoOffset_s;
[peakAmplitude, linearIdx] = max(abs(audioSegment(:)));
[peakIdx, peakChannel] = ind2sub(size(audioSegment), linearIdx);
peakTime = timeVector(peakIdx);
fprintf('Loudest sample (not confirmed onset): %.6f video seconds, channel %d, |amplitude| %.4f\n', ...
    peakTime, peakChannel, peakAmplitude);

figure;
layout = tiledlayout(3, 1);
axWave = nexttile;
plot(timeVector, audioSegment);
title(sprintf('%s: audio inspection (all %d channels)', trialName, numChannels));
ylabel('Signed amplitude');
grid on;
hold on;
plot(peakTime, audioSegment(peakIdx, peakChannel), 'ro', 'MarkerSize', 8);
xline(peakTime, '--r', sprintf('Loudest: %.4f s', peakTime));
hold off;

% Short RMS bins reveal sustained sound changes without phase cancellation
% between channels. This is a display aid, not an automatic onset threshold.
axRms = nexttile;
binSamples = max(1, round(0.02 * fs));
binStarts = (1:binSamples:size(audioSegment, 1))';
rmsData = zeros(numel(binStarts), numChannels);
rmsTime = zeros(numel(binStarts), 1);
for bin = 1:numel(binStarts)
    idx = binStarts(bin):min(binStarts(bin) + binSamples - 1, size(audioSegment, 1));
    rmsData(bin, :) = sqrt(mean(audioSegment(idx, :).^2, 1));
    rmsTime(bin) = mean(timeVector(idx));
end
plot(rmsTime, rmsData);
ylabel('20 ms RMS');
legend(compose('Channel %d', 1:numChannels), 'Location', 'best');
grid on;

axAngle = nexttile;
angleFile = fullfile(repoDir, 'Data_Generation', 'data', [trialName '_data.csv']);
if isfile(angleFile)
    angleData = readmatrix(angleFile);
    plot(angleData(:, 1), angleData(:, 2));
    ylabel('Angle (degrees)');
else
    title('No angular-position CSV found');
end
grid on;
xlabel('Original video time (s)');
linkaxes([axWave axRms axAngle], 'x');
xlim(axWave, [timeVector(1), endTime + audioToVideoOffset_s]);

if playbackEnabled
    validateattributes(playbackChannels, {'numeric'}, ...
        {'vector', 'integer', 'positive', '<=', numChannels, 'nonempty'});
    assert(numel(playbackChannels) <= 2, 'Choose at most two channels for playback.');
    validateattributes(playbackGain, {'numeric'}, {'scalar', 'real', 'finite', '>=', 0, '<=', 1});
    player = audioplayer(playbackGain * audioSegment(:, playbackChannels), fs);
    hLine = xline(axWave, timeVector(1), '-r');
    disp('Playing selected channels at reduced gain.');
    play(player);
    while isplaying(player) && isgraphics(hLine)
        hLine.Value = timeVector(1) + (player.CurrentSample - 1) / fs;
        drawnow limitrate;
        pause(0.02);
    end
    stop(player);
end

% Use confirmed sound cues to help choose integration windows in main.py.
% startup_window / shutdown_window are TORQUE integration bounds, not audio
% search bounds. Choose them on the reconstructed torque and check that its
% local cumulative integral settles as the bounds are varied.
