function candidates = estimateSoundTimes(videoFileName, settings, angleData)
%ESTIMATESOUNDTIMES Estimate pump-event cues using the method from the former sound.m script.
%   candidates = estimateSoundTimes(videoFileName, config.audio) returns a
%   struct with startup and shutdown times in original VIDEO seconds.
%   Each search interval is supplied in AUDIO seconds; the configured offset
%   converts a sample timestamp using video time = audio time + offset.
%
%   The loudest sample is a candidate cue, not a confirmed onset. Generation
%   saves these candidates; manual corrections are kept separately in the
%   trial metadata. This helper does not read or overwrite those corrections.
%   Inspection and playback are optional and do not alter the calculation.
%   Optional angleData is [original video seconds, angle degrees] and restores
%   the original angular-position panel alongside the audio waveform and RMS.
    if nargin < 3
        angleData = [];
    end
    % Read all channels once. Playback channel choices must not change which
    % sample is selected for either event.
    [audioData, fs] = audioread(videoFileName);
    totalSamples = size(audioData, 1);
    numChannels = size(audioData, 2);
    events = {'startup', 'shutdown'};
    for eventIndex = 1:numel(events)
        event = events{eventIndex};
        % Startup and shutdown have independent searches; shutdown is not
        % inferred by adding a fixed pump duration to the startup candidate.
        interval = settings.([event '_search_s']);
        startTime = interval(1);
        endTime = min(interval(2), totalSamples / fs);
        validateattributes(startTime, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
        validateattributes(endTime, {'numeric'}, {'scalar', 'real', 'finite', '>', startTime});
        % Half-open interval [startTime, endTime), as in the original calculation.
        startIndex = ceil(startTime * fs) + 1;
        endIndex = min(totalSamples, ceil(endTime * fs));
        assert(endIndex >= startIndex, '%s search contains no audio samples.', event);
        segment = audioData(startIndex:endIndex, :);
        time = ((startIndex:endIndex)' - 1) / fs + settings.audio_to_video_offset_s;
        % Flattening searches every sample/channel pair without averaging
        % channels, which could cancel an otherwise strong sound cue.
        [amplitude, linearIndex] = max(abs(segment(:)));
        [peakIndex, peakChannel] = ind2sub(size(segment), linearIndex);
        candidates.(event) = time(peakIndex);
        fprintf('%s loudest sample (candidate): %.6f video s, channel %d, amplitude %.4f\n', ...
            event, time(peakIndex), peakChannel, amplitude);
        if settings.inspect || settings.playback_enabled
            % All panels use video seconds so they can be compared with the
            % angle CSV and the integration-window selections in Python.
            figure;
            tiledlayout(3, 1);
            axWave = nexttile;
            plot(time, segment);
            hold on;
            plot(time(peakIndex), segment(peakIndex, peakChannel), 'ro', 'MarkerSize', 8);
            xline(time(peakIndex), '--r', sprintf('Loudest: %.4f s', time(peakIndex)));
            hold off;
            title(sprintf('%s: %s candidate', videoFileName, event), 'Interpreter', 'none');
            ylabel('Signed amplitude');
            grid on;
            % Short RMS bins reveal sustained sound changes without phase cancellation
            % between channels. This is a display aid, not an automatic onset threshold.
            axRms = nexttile;
            binSamples = max(1, round(0.02 * fs));
            binStarts = (1:binSamples:size(segment, 1))';
            rmsData = zeros(numel(binStarts), numChannels);
            rmsTime = zeros(numel(binStarts), 1);
            for bin = 1:numel(binStarts)
                idx = binStarts(bin):min(binStarts(bin) + binSamples - 1, size(segment, 1));
                rmsData(bin, :) = sqrt(mean(segment(idx, :).^2, 1));
                rmsTime(bin) = mean(time(idx));
            end
            plot(rmsTime, rmsData);
            ylabel('20 ms RMS');
            legend(compose('Channel %d', 1:numChannels), 'Location', 'best');
            grid on;

            axAngle = nexttile;
            if ~isempty(angleData)
                plot(angleData(:, 1), angleData(:, 2));
                ylabel('Angle (degrees)');
            else
                title('No angular-position data supplied');
            end
            grid on;
            xlabel('Original video time (s)');
            linkaxes([axWave axRms axAngle], 'x');
            xlim(axWave, [time(1), endTime + settings.audio_to_video_offset_s]);

            if settings.playback_enabled
                % A fixed gain attenuates playback without normalizing it.
                % Retain the original moving cursor while each excerpt plays.
                channels = settings.playback_channels;
                validateattributes(channels, {'numeric'}, ...
                    {'vector', 'integer', 'positive', '<=', numChannels, 'nonempty'});
                assert(numel(channels) <= 2, 'Choose at most two playback channels.');
                validateattributes(settings.playback_gain, {'numeric'}, ...
                    {'scalar', 'real', 'finite', '>=', 0, '<=', 1});
                player = audioplayer(settings.playback_gain * segment(:, channels), fs);
                hLine = xline(axWave, time(1), '-r');
                disp('Playing selected channels at reduced gain.');
                play(player);
                while isplaying(player) && isgraphics(hLine)
                    hLine.Value = time(1) + (player.CurrentSample - 1) / fs;
                    drawnow limitrate;
                    pause(0.02);
                end
                stop(player);
            end
        end
    end
end
