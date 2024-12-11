function openephysDig = loadOpenEphysDigitalNidaq(session, varargin)
% function to load digital inputs from OpenEphys NI-DAQ Board and align the data properly
%
% INPUT
% session: session struct
% Optional inputs:
%   channelNum: Channel number for NI-DAQ board (0 for clock, 1 for camera)
%   probeLetter: 'A' or 'B' for specific probe timestamps. If empty, uses ProbeA
%
% OUTPUTS
% openephysDig.on: on-state changes for each channel [in seconds]
% openephysDig.off: off-state changes for each channel [in seconds]
%
% By Mingze Dou
% mingzedou.mail@gmail.com

p = inputParser;
addParameter(p,'channelNum', [], @isnumeric);
addParameter(p,'probeLetter', '', @ischar);
parse(p,varargin{:});
parameters = p.Results;

TTL_paths = {};
epochs_startTime = [];
ephys_t0 = [];

% Determine paths and timestamps for each epoch
for i = 1:numel(session.epochs)
    TTL_paths{i} = fullfile(session.epochs{i}.name, 'events', 'NI-DAQmx-105.PXIe-6341', 'TTL');
    disp(['Checking TTL path for epoch ' num2str(i) ': ' fullfile(session.general.basePath, TTL_paths{i})])
    
    % Set timestamp path based on probe selection
    if ~isempty(parameters.probeLetter)
        probeStr = parameters.probeLetter;
    else
        probeStr = 'A';
    end
    
    timestampPath = fullfile(session.general.basePath, session.epochs{i}.name, 'continuous', ...
        ['Neuropix-PXI-100.Probe' probeStr], 'timestamps.npy');
    disp(['Checking timestamp path for epoch ' num2str(i) ': ' timestampPath])
    
    if ~exist(timestampPath, 'file')
        error(['Timestamp file not found for Probe' probeStr ' in epoch ' num2str(i)]);
    end
    
    temp = readNPY(timestampPath);
    disp(['Found ' num2str(length(temp)) ' timestamps, first value: ' num2str(temp(1))])
    
    epochs_startTime(i) = session.epochs{i}.startTime;
    ephys_t0(i) = double(temp(1));
    disp(['ephys_t0 for epoch ' num2str(i) ': ' num2str(ephys_t0(i))])
end

% Initialize output structure
openephysDig = {};

% Load and process first epoch data
basePath = fullfile(session.general.basePath, TTL_paths{1});

% Load raw data
timestamps = readNPY(fullfile(basePath,'timestamps.npy'));
disp(['Found ' num2str(length(timestamps)) ' TTL timestamps'])
full_words = readNPY(fullfile(basePath,'full_words.npy'));
channel_states = readNPY(fullfile(basePath,'states.npy'));
channels = readNPY(fullfile(basePath,'sample_numbers.npy'));

% Filter by channel first if specified
if ~isempty(parameters.channelNum)
    channel_mask = bitand(full_words, 2^parameters.channelNum) > 0;
    timestamps = timestamps(channel_mask);
    channel_states = channel_states(channel_mask);
    channels = channels(channel_mask);
    full_words = full_words(channel_mask);
end

% Then store in structure with alignment
openephysDig.timestamps = epochs_startTime(1) + double(timestamps) - ephys_t0(1);
openephysDig.channel_states = channel_states;
openephysDig.channels = channels;
openephysDig.full_words = full_words;

openephysDig.on{1} = double(openephysDig.timestamps(openephysDig.channel_states == 1));
openephysDig.off{1} = double(openephysDig.timestamps(openephysDig.channel_states == -1));

% Process additional epochs if present
if length(TTL_paths) > 1
    openephysDig.nTimestampsPrFile(1) = numel(openephysDig.timestamps);
    openephysDig.nOnPrFile(1) = numel(openephysDig.on{1});
    openephysDig.nOffPrFile(1) = numel(openephysDig.off{1});
    
    for i = 2:length(TTL_paths)
        basePath = fullfile(session.general.basePath, TTL_paths{i});

        timestamps = readNPY(fullfile(basePath,'timestamps.npy'));
        full_words = readNPY(fullfile(basePath,'full_words.npy'));
        channel_states = readNPY(fullfile(basePath,'states.npy'));
        channels = readNPY(fullfile(basePath,'sample_numbers.npy'));
        
        % Filter for specific channel if requested
        if ~isempty(parameters.channelNum)
            channel_mask = bitand(full_words, 2^parameters.channelNum) > 0;
            timestamps = timestamps(channel_mask);
            channel_states = channel_states(channel_mask);
            channels = channels(channel_mask);
            full_words = full_words(channel_mask);
        end

        timestamps = epochs_startTime(i) + double(timestamps) - ephys_t0(i);
        
        openephysDig.timestamps = [openephysDig.timestamps; timestamps];
        openephysDig.channel_states = [openephysDig.channel_states; channel_states];
        openephysDig.channels = [openephysDig.channels; channels];
        openephysDig.full_words = [openephysDig.full_words; full_words];
        
        openephysDig.on{1} = [openephysDig.on{1}; double(timestamps(channel_states == 1))];
        openephysDig.off{1} = [openephysDig.off{1}; double(timestamps(channel_states == -1))];
        openephysDig.nTimestampsPrFile(i) = numel(timestamps);
        openephysDig.nOnPrFile(i) = sum(channel_states == 1);
        openephysDig.nOffPrFile(i) = sum(channel_states == -1);
    end
end

% Attaching info about how the data was processed
openephysDig.processinginfo.function = 'loadOpenEphysDigitalNidaq';
openephysDig.processinginfo.version = 1;
openephysDig.processinginfo.date = datetime('now');
openephysDig.processinginfo.params.basepath = session.general.basePath;
openephysDig.processinginfo.params.basename = session.general.name;
openephysDig.processinginfo.params.TTL_paths = TTL_paths;

try
    openephysDig.processinginfo.username = char(java.lang.System.getProperty('user.name'));
    openephysDig.processinginfo.hostname = char(java.net.InetAddress.getLocalHost.getHostName);
catch
    disp('Failed to retrieve system info.')
end

% Saving data
saveStruct(openephysDig,'digitalseries','session',session);
end