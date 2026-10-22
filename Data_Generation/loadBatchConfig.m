function config = loadBatchConfig(configFile)
%LOADBATCHCONFIG Read the batch settings shared by MATLAB and Python.
%   config = loadBatchConfig() reads experiment_config.json at the repo root.
%   config = loadBatchConfig(configFile) reads a specified JSON configuration.
%   Returns the decoded settings, with config.paths converted to full paths
%   when configFile is absolute, and config.trials expanded to full trial IDs.
%
%   List trial_groups such as ["1000f", "500r"] and set trials_per_group once.
%   A count of 3 produces 1000f1..1000f3, then 500r1..500r3, in that order.
%   An array [2, 3] instead gives 2 repeats of 1000f and 3 repeats of 500r.
%   Array counts must have one positive integer for each group, in group order.
%   The expanded config.trials exists only in memory; no trial manifest is saved.
%   Edit the JSON to change groups, repeat count, paths, video geometry,
%   laser thresholds, physical constants, audio searches or analysis settings.
%   One config describes one compatible camera/ruler batch; spring calibration
%   keeps its separate setup values in calibrateSpring.m.
    if nargin < 1
        % Find the repo from this helper, so callers can run from another folder.
        repoDir = fileparts(fileparts(mfilename('fullpath')));
        configFile = fullfile(repoDir, 'experiment_config.json');
    end
    configText = fileread(configFile);
    config = jsondecode(configText);
    configDir = fileparts(configFile);
    names = fieldnames(config.paths);
    % Match Python's load_config: relative paths belong to the JSON's folder.
    % Leave Unix absolute paths and Windows drive-letter paths unchanged.
    for i = 1:numel(names)
        value = string(config.paths.(names{i}));
        if ~startsWith(value, filesep) && isempty(regexp(value, '^[A-Za-z]:[\\/]', 'once'))
            value = fullfile(configDir, value);
        end
        config.paths.(names{i}) = value;
    end
    % JSON string arrays decode as cells. Require {Re}{f/r} group names,
    % without a repeat number. Reject malformed settings before any
    % script constructs filenames or removes old generated outputs.
    assert(iscellstr(config.trial_groups) && ~isempty(config.trial_groups), ...
        'trial_groups must be a nonempty array of group IDs such as 1000f.');
    config.trial_groups = string(config.trial_groups);
    assert(numel(unique(config.trial_groups)) == numel(config.trial_groups), ...
        'Configured trial groups must be unique.');
    for group = config.trial_groups(:)'
        assert(~isempty(regexp(group, '^[1-9]\d*[fr]$', 'once')), ...
            'Invalid trial group: %s; expected e.g. 1000f.', group);
    end
    counts = config.trials_per_group;
    validateattributes(counts, {'numeric'}, ...
        {'vector', 'nonempty', 'real', 'finite', 'integer', 'positive'}, mfilename, 'trials_per_group');
    % jsondecode turns both 1 and [1] into scalar doubles. Preserve the JSON
    % distinction so an array still requires exactly one count per group.
    if isempty(regexp(configText, '"trials_per_group"\s*:\s*\[', 'once'))
        assert(isscalar(counts), 'trials_per_group must be a positive integer or an array.');
        counts = repmat(counts, numel(config.trial_groups), 1);
    else
        % MATLAB also flattens nested singleton arrays; reject those explicitly
        % so only a flat list of counts is accepted, matching Python.
        assert(isempty(regexp(configText, '"trials_per_group"\s*:\s*\[[^\]]*\[', 'once')), ...
            'trials_per_group must be a flat array of positive integers.');
        assert(numel(counts) == numel(config.trial_groups), ...
            'trials_per_group must contain one positive integer per trial group.');
    end

    % Match Python's load_config: expand all repeats of the first group,
    % followed by all repeats of the next. Existing scripts still loop over
    % config.trials, and the first group's repeat 1 supplies the batch reference.
    % No folder discovery is done; missing recordings are reported by the caller.
    config.trials = strings(sum(counts), 1);
    nextTrial = 1;
    for groupIndex = 1:numel(config.trial_groups)
        group = config.trial_groups(groupIndex);
        for trialNumber = 1:counts(groupIndex)
            config.trials(nextTrial) = group + string(trialNumber);
            nextTrial = nextTrial + 1;
        end
    end
end
