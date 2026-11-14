function [config, trialPath] = loadBatchConfig(configFile)
%LOADBATCHCONFIG Shared paths and inclusive ingestion ranges (never analysis IDs).
    if nargin < 1
        configFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'experiment_config.json');
    end
    configText = fileread(configFile);
    config = jsondecode(configText);
    configDir = fileparts(configFile);
    if isempty(configDir), configDir = pwd; end
    names = fieldnames(config.paths);
    for i = 1:numel(names)
        value = string(config.paths.(names{i}));
        if ~startsWith(value, filesep) && isempty(regexp(value, '^[A-Za-z]:[\\/]', 'once'))
            value = fullfile(configDir, value);
        end
        config.paths.(names{i}) = value;
    end
    % Read original JSON keys: jsondecode makes "1000r" into "x1000r".
    % Strict flat grammar also rejects nested arrays, scalars and duplicate keys.
    block = regexp(configText, '"trial_ranges"\s*:\s*\{([^{}]*)\}', 'tokens');
    assert(numel(block) == 1, 'trial_ranges must be one object.');
    body = strjoin(block{1}, '');
    entry = '"([1-9]\d*[fr])"\s*:\s*(\[[^\[\]]*\])';
    grammar = ['^\s*(' entry '(\s*,\s*' entry ')*)?\s*$'];
    assert(isempty(strtrim(body)) || ~isempty(regexp(body, grammar, 'once')), 'Invalid trial_ranges group or range.');
    entries = regexp(body, entry, 'tokens');
    groups = strings(0,1);
    config.ingestion_trials = strings(0,1);
    for i = 1:numel(entries)
        group = string(entries{i}{1});
        assert(~any(groups == group), 'Duplicate trial_ranges group: %s', group);
        groups(end+1,1) = group;
        raw = entries{i}{2};
        if ~isempty(regexp(raw, '^\[\s*\]$', 'once')), continue; end
        % Exactly two JSON numeric tokens; booleans and nested singleton arrays
        % must not become numeric endpoints through MATLAB's decoding rules.
        number = '-?(0|[1-9]\d*)(\.\d+)?([eE][+-]?\d+)?';
        assert(~isempty(regexp(raw, ['^\[\s*' number '\s*,\s*' number '\s*\]$'], 'once')), ...
            'Ranges require [] or [first,last].');
        bounds = jsondecode(raw);
        validateattributes(bounds, {'double'}, {'numel',2,'real','finite','integer','positive','<=',flintmax-1});
        assert(bounds(1) <= bounds(2), 'Reversed trial range: %s', group);
        config.ingestion_trials = [config.ingestion_trials; group + string((bounds(1):bounds(2))')]; %#ok<AGROW>
    end
    if ~isfield(config.video, 'regenerate_trials'), config.video.regenerate_trials = {}; end
    ids = config.video.regenerate_trials;
    assert((isempty(ids) && isnumeric(ids)) || iscellstr(ids), 'regenerate_trials must be an array of trial IDs.');
    for i = 1:numel(ids)
        assert(any(config.ingestion_trials == string(ids{i})), 'Regeneration IDs must be in ingestion ranges.');
    end
    config.video.regenerate_trials = string(ids);
    trialPath = @storagePath;
end

function path = storagePath(root, trial, suffix)
% Grouped destination with read compatibility for existing flat artifacts.
    group = regexp(char(trial), '^[1-9]\d*[fr](?=[1-9]\d*$)', 'match', 'once');
    assert(~isempty(group), 'Invalid trial ID: %s', trial);
    grouped = fullfile(root, group, string(trial) + suffix);
    legacy = fullfile(root, string(trial) + suffix);
    assert(~(isfile(grouped) || isfolder(grouped)) || ~(isfile(legacy) || isfolder(legacy)), ...
        'Trial storage collision: %s and %s', grouped, legacy);
    if isfile(legacy) || isfolder(legacy), path = legacy; else, path = grouped; end
end
