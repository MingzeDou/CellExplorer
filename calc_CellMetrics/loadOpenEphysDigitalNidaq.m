function openephysDig = loadOpenEphysDigitalNidaq(session, varargin)
p = inputParser;
addParameter(p,'channelNum', [], @isnumeric);
addParameter(p,'probeLetter', 'A', @ischar);
parse(p,varargin{:});
parameters = p.Results;

% Initialize
ttlPaths = {};
epochsStartTime = [];
ephysT0 = [];
validEpochs = [];
openephysDig.on = {[]};
openephysDig.off = {[]};
openephysDig.timestamps = [];
openephysDig.states = [];
openephysDig.epochNum = [];
openephysDig.diagnostics = struct();

% Find valid epochs
for i = 1:numel(session.epochs)
    ttlPath = fullfile(session.epochs{i}.name, 'events', 'NI-DAQmx-116.PXIe-6341', 'TTL');
    fullWordsPath = fullfile(session.general.basePath, ttlPath, 'full_words.npy');
    
    if exist(fullWordsPath, 'file') && (~isempty(parameters.channelNum) && ...
            any(bitand(readNPY(fullWordsPath), 2^parameters.channelNum)) || isempty(parameters.channelNum))
        validEpochs = [validEpochs i];
        ttlPaths{end+1} = ttlPath;
        
        % Get probe timestamps
        probeStr = parameters.probeLetter;
        if ~isempty(parameters.probeLetter), probeStr = parameters.probeLetter; end
        
        timestampPath = fullfile(session.general.basePath, session.epochs{i}.name, 'continuous', ...
            ['Neuropix-PXI-103.Probe' probeStr], 'timestamps.npy');
        
        if ~exist(timestampPath, 'file')
            error(['Timestamp file not found for Probe' probeStr ' in epoch ' num2str(i)]);
        end
        
        timestampData = readNPY(timestampPath);
        epochsStartTime(end+1) = session.epochs{i}.startTime;
        ephysT0(end+1) = double(timestampData(1));
    end
end

if isempty(validEpochs)
    error('No valid epochs found with the specified channel');
end

% Process epochs
for i = 1:length(validEpochs)
    currentEpoch = validEpochs(i);
    basePath = fullfile(session.general.basePath, ttlPaths{i});
    timestamps = readNPY(fullfile(basePath,'timestamps.npy'));
    fullWords = readNPY(fullfile(basePath,'full_words.npy'));
    if isempty(timestamps) || isempty(fullWords)
        warning('Empty data found in epoch %d', currentEpoch);
        continue;
    end

    % Align timestamps and get states
    alignedTimestamps = epochsStartTime(i) + double(timestamps) - ephysT0(i);
    channel_states = bitget(fullWords, parameters.channelNum + 1);
    [uniqueTimestamps, uniqueIdx] = unique(alignedTimestamps);
    uniqueStates = channel_states(uniqueIdx);

    % Initialize arrays for this epoch
    rising_edges = [];
    falling_edges = [];

    % Process the state sequence
    state_idx = 1;
    while state_idx < length(uniqueStates)
        % Find next rising edge
        while state_idx < length(uniqueStates) && uniqueStates(state_idx) == 0
            state_idx = state_idx + 1;
        end

        if state_idx < length(uniqueStates)
            % Found a rising edge
            rising_edge_idx = state_idx;

            % Look for the next falling edge for this specific channel
            state_idx = state_idx + 1;
            while state_idx < length(uniqueStates) && uniqueStates(state_idx) == 1
                state_idx = state_idx + 1;
            end

            if state_idx < length(uniqueStates)
                % Found the corresponding falling edge
                falling_edge_idx = state_idx;

                % Store the edges
                rising_edges = [rising_edges; rising_edge_idx];
                falling_edges = [falling_edges; falling_edge_idx];
            else
                warning('Found rising edge without corresponding falling edge at timestamp %f in epoch %d', ...
                    uniqueTimestamps(rising_edge_idx), currentEpoch);
            end
        end
    end

    % Store results using the correct timestamps
    openephysDig.on{1} = [openephysDig.on{1}; uniqueTimestamps(rising_edges)];
    openephysDig.off{1} = [openephysDig.off{1}; uniqueTimestamps(falling_edges)];
    openephysDig.timestamps = [openephysDig.timestamps; uniqueTimestamps];
    openephysDig.states = [openephysDig.states; uniqueStates];
    openephysDig.epochNum = [openephysDig.epochNum; repmat(currentEpoch, length(uniqueTimestamps), 1)];
end
% Save
saveStruct(openephysDig,'digitalseries','session',session);
end