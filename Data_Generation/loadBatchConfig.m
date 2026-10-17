function config = loadBatchConfig(configFile)
%LOADBATCHCONFIG Read the batch settings shared by MATLAB and Python.
%   config = loadBatchConfig() reads experiment_config.json at the repo root.
%   config = loadBatchConfig(configFile) reads a specified JSON configuration.
%   Returns the decoded settings, with config.paths converted to full paths
%   when configFile is absolute, and config.trials expanded to full trial IDs.
%
%   List trial_groups such as ["1000f", "500r"] and set trials_per_group once.
%   A count of 3 produces 1000f1..1000f3, then 500r1..500r3, in that order.
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
    config = jsondecode(fileread(configFile));
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
    validateattributes(config.trials_per_group, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'positive'}, mfilename, 'trials_per_group');

    % Match Python's load_config: expand all repeats of the first group,
    % followed by all repeats of the next. Existing scripts still loop over
    % config.trials, and the first group's repeat 1 supplies the batch reference.
    % No folder discovery is done; missing recordings are reported by the caller.
    config.trials = strings(numel(config.trial_groups) * config.trials_per_group, 1);
    nextTrial = 1;
    for group = config.trial_groups(:)'
        for trialNumber = 1:config.trials_per_group
            config.trials(nextTrial) = group + string(trialNumber);
            nextTrial = nextTrial + 1;
        end
    end
end
