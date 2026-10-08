% force using the default audio output
% Please don't wear headphones on!!!
info = audiodevinfo;
disp('Default audio output device:');
disp(info.output(1));

% read audio file
repoDir = fileparts(mfilename('fullpath'));
trialName = '1000r1';
videoFileName = fullfile(repoDir, 'Data_Generation', 'trials', [trialName '.MTS']);
disp(['Reading file: ', videoFileName, ' ...']);
[audioData, fs] = audioread(videoFileName);

% check channels
numChannels = size(audioData, 2);
disp(['Number of audio channels: ', num2str(numChannels)]);

% choose the time range you want 
startTime = 9.5; 
endTime = 13;   

% sample indices
startIndex = floor(startTime * fs) + 1;
endIndex = floor(endTime * fs);
totalSamples = length(audioData); 
if endIndex > totalSamples
    endIndex = totalSamples;
end

% extract segment
% the reason why i only take at most 2 tracks is because of the later
% function audioplayer take no more than 2 tracks
audioSegment = audioData(startIndex:endIndex, :);
if numChannels > 2
    disp('Using only the first two channels for playback');
    audioSegment = audioSegment(:, 1:2); 
end

% Compute maximum absolute amplitude across channels
if size(audioSegment, 2) == 2
    absData = max(abs(audioSegment), [], 2); 
else
    absData = abs(audioSegment);
end
[peakAmplitude, peakIdx] = max(absData);     
peakTime = startTime + (peakIdx - 1) / fs;    % convert to seconds

% print peak time information
fprintf('A possible starting time according to the audio is %.2f seconds\n', peakTime);


% plot waveform
timeVector = linspace(startTime, endTime, length(audioSegment));
figure; 
plot(timeVector, audioSegment);
title(['Audio Waveform: ', num2str(startTime), 's to ', num2str(endTime), 's']);
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;
hold on;
plot(peakTime, peakAmplitude, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xline(peakTime, '--r', ['Peak: ', num2str(peakTime, '%.2f'), 's']);
if size(audioSegment, 2) == 2
    legend('Channel 1', 'Channel 2', 'Peak Position');
else
    legend('Channel 1', 'Peak Position');
end

% progress line
yl = ylim;
hLine = plot([startTime startTime], yl, 'r', 'LineWidth', 2);
xlabel('Time (s)');

% playback
player = audioplayer(audioSegment, fs);
disp(['Playing from ', num2str(startTime), 's to ', num2str(endTime), 's ...']);
play(player);

% animation loop
fprintf('Progress: ');
while isplaying(player)
    currentSec = startTime + (player.CurrentSample / fs);
    set(hLine, 'XData', [currentSec currentSec]);
    progress = (currentSec - startTime) / (endTime - startTime) * 100;
    fprintf('\b\b\b\b\b\b%5.1f%%', progress);
    drawnow limitrate;
end
hold off;