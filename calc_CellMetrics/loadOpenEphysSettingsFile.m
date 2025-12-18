function session = loadOpenEphysSettingsFile(file1,session,varargin)
% Loading structure.oebin -  a json structure created by OpenEphys
% https://open-ephys.github.io/gui-docs/User-Manual/Recording-data/Binary-format.html

p = inputParser;
addParameter(p,'probeLetter','A',@ischar);
parse(p,varargin{:})
parameters = p.Results;

switch parameters.probeLetter
    case 'A'
        stream_id = 1;
    case 'B'
        stream_id = 2;
    case 'C'
        stream_id = 3;
    otherwise
        stream_id = 1;
end

if ~exist(file1,'file')
    warning(['Failed to load Open Ephys settings: ', file1])
    return
end

disp(['Loading Open Ephys settings: ', file1])
text = fileread(file1);
openEphys_metadata = jsondecode(text);

% Handle "GUIVersion" vs "GUI version"
if     isfield(openEphys_metadata,'GUIVersion')
    GUIVersion = openEphys_metadata.GUIVersion;
elseif isfield(openEphys_metadata,'GUI version')
    GUIVersion = openEphys_metadata.('GUI version');
else
    GUIVersion = '';
end

% ===================== OE 0.6.7 (unchanged) =====================
if strcmp(GUIVersion, '0.6.7')
    % Importing metadata
    session.extracellular.sr = openEphys_metadata.continuous(stream_id).sample_rate;
    session.extracellular.nChannels = openEphys_metadata.continuous(stream_id).num_channels;
    session.extracellular.equipment = 'OpenEpys Neuropix-PXI';
    session.extracellular.leastSignificantBit = 0.195;
    session.extracellular.fileFormat = 'dat';
    session.extracellular.precision = 'int16';
    session.extracellular.fileName = '';
    
    % Electrode groups and channel mapping
    channelmapping = [];
    for i = 1:session.extracellular.nChannels
        if isstruct(openEphys_metadata.continuous(stream_id).channels) 
            if isfield(openEphys_metadata.continuous(stream_id).channels(i),'channel_metadata')
                channelmapping(i) = openEphys_metadata.continuous(stream_id).channels(i).channel_metadata.value+1;
            end
        elseif isfield(openEphys_metadata.continuous(stream_id).channels{i},'channel_metadata')
            channelmapping(i) = openEphys_metadata.continuous(stream_id).channels{i}.channel_metadata.value+1;
        end
    end

    chanCoords = generateChanCoords_Neuropixel2;

    session.extracellular.chanCoords = chanCoords;
    session.extracellular.chanCoords.x = chanCoords.x(channelmapping);
    session.extracellular.chanCoords.y = chanCoords.y(channelmapping);
    session.extracellular.chanCoords.shank = chanCoords.shank(channelmapping);
    session.extracellular.chanCoords.channels = chanCoords.channels(channelmapping);

    session.extracellular.electrodeGroups.channels = {};
    
    shanks = unique(session.extracellular.chanCoords.shank);
    for j = 1:length(shanks)
        session.extracellular.electrodeGroups.channels{j} = find(session.extracellular.chanCoords.shank == shanks(j));
        [~,idx] = sortrows([...
            session.extracellular.chanCoords.x(session.extracellular.electrodeGroups.channels{j});...
            session.extracellular.chanCoords.y(session.extracellular.electrodeGroups.channels{j})...
            ]',[2,1],{'descend','ascend'});
        session.extracellular.electrodeGroups.channels{j} = session.extracellular.electrodeGroups.channels{j}(idx);
        session.extracellular.electrodeGroups.label{j} = ['shanks',num2str(shanks(j))];
    end
    session.extracellular.nElectrodeGroups = length(shanks);

    figure
    plot(chanCoords.x,chanCoords.y,'.k'), hold on
    site_cmap = hot(session.extracellular.nElectrodeGroups);
    for j = 1:session.extracellular.nElectrodeGroups
        x = session.extracellular.chanCoords.x(session.extracellular.electrodeGroups.channels{j});
        y = session.extracellular.chanCoords.y(session.extracellular.electrodeGroups.channels{j});
        plot(x,y,'s',MarkerFaceColor=site_cmap(j,:))
    end
    xlabel('X position (µm)'), ylabel('Y position (µm)'), title('Neuropixel site selection', file1)
    axis equal

    session.extracellular.spikeGroups.channels = session.extracellular.electrodeGroups.channels;
    session.extracellular.nSpikeGroups = session.extracellular.nElectrodeGroups;

    % LFP data stream
    if length(openEphys_metadata.continuous)>1
        session.extracellular.srLfp = openEphys_metadata.continuous(2).sample_rate;
    end

    % Loading events from timestamps
    session.timeSeries.dig.fileName = [openEphys_metadata.events{1}.folder_name,'timestamps.npy'];
    session.timeSeries.dig.fileFormat = 'npy';
    session.timeSeries.dig.precision  = openEphys_metadata.events{1}.type;
    session.timeSeries.dig.nChannels = 1;
    session.timeSeries.dig.sr = openEphys_metadata.events{1}.sample_rate;
    session.timeSeries.dig.equipment = 'OpenEpys Neuropix-PXI';
    session.timeSeries.dig.nSamples = [];
    session.timeSeries.dig.leastSignificantBit = [];
end

% ===================== OE 1.0.1 (NP2.0 or NP1.0) =====================
if strcmp(GUIVersion, '1.0.1')
    % Importing metadata (unchanged header)
    session.extracellular.sr = openEphys_metadata.continuous(stream_id).sample_rate;
    session.extracellular.nChannels = openEphys_metadata.continuous(stream_id).num_channels;
    session.extracellular.equipment = 'Open Ephys Neuropix-PXI';
    session.extracellular.leastSignificantBit = 0.195;
    session.extracellular.fileFormat = 'dat';
    session.extracellular.precision = 'int16';
    session.extracellular.fileName = '';

    % --- Extract electrode_index for detection/mapping ---
    cstr = openEphys_metadata.continuous(stream_id);
    chList = cstr.channels;
    if isstruct(chList)
        n = numel(chList);
    elseif iscell(chList)
        n = numel(chList);
    else
        error('Unexpected type for continuous(stream_id).channels');
    end
    if n ~= session.extracellular.nChannels
        warning('Metadata nChannels (%d) != channels entries (%d). Using channels entries.', ...
            session.extracellular.nChannels, n);
        session.extracellular.nChannels = n;
    end

    electrode_index_vec = nan(1,n);
    for i = 1:n
        if isstruct(chList), c = chList(i); else, c = chList{i}; end
        ei = [];
        if isfield(c,'channel_metadata') && ~isempty(c.channel_metadata)
            meta = c.channel_metadata;
            if iscell(meta), meta = [meta{:}]; end
            if isstruct(meta)
                if numel(meta)==1 && isfield(meta,'value') && ~isempty(meta.value)
                    v = meta.value;
                    if isnumeric(v), ei = double(v(1)); end
                else
                    k = find(arrayfun(@(m) isfield(m,'name') && strcmpi(m.name,'electrode_index'), meta), 1);
                    if ~isempty(k) && isfield(meta(k),'value') && ~isempty(meta(k).value)
                        v = meta(k).value;
                        if isnumeric(v), ei = double(v(1)); end
                    end
                end
            end
        end
        if isempty(ei) && isfield(c,'metadata') && isfield(c.metadata,'electrodeID')
            ei = double(c.metadata.electrodeID);
        end
        if isempty(ei)
            ei = i-1; % fallback assume 0-indexed acquisition order
        end
        electrode_index_vec(i) = ei;
    end

    % --- Detect NP1.0 vs NP2.0 by electrode_index range ---
    ei_max = max(electrode_index_vec);
    isNP10 = (ei_max <= 959);  % NP1.0 ≤ 960 sites; else NP2.0

    if isNP10
        % ================= NP1.0 path =================
        chanCoords = generateChanCoords_Neuropixel1;
        geomLen = numel(chanCoords.x);

        % Normalize indexing to MATLAB 1-based
        if any(electrode_index_vec == 0) || max(electrode_index_vec) == geomLen-1
            ei1 = electrode_index_vec + 1;
        elseif max(electrode_index_vec) <= geomLen
            ei1 = electrode_index_vec;
        else
            error('electrode_index out of geometry bounds: max=%d > %d', max(electrode_index_vec), geomLen);
        end

        % Bottom→top physical order by electrode_index
        [~, order] = sort(ei1, 'ascend');  % bottom-up

        % Apply geometry by physical index
        session.extracellular.chanCoords = chanCoords;
        session.extracellular.chanCoords.x        = chanCoords.x(ei1);
        session.extracellular.chanCoords.y        = chanCoords.y(ei1);
        session.extracellular.chanCoords.shank    = ones(1, numel(ei1));      % single shank
        session.extracellular.chanCoords.channels = chanCoords.channels(ei1);

        % Electrode/spike groups: acquisition channels arranged bottom-up
        session.extracellular.electrodeGroups.channels = {order(:).'};
        session.extracellular.electrodeGroups.label    = {'shank1'};
        session.extracellular.nElectrodeGroups         = 1;

        session.extracellular.spikeGroups.channels     = session.extracellular.electrodeGroups.channels;
        session.extracellular.nSpikeGroups             = 1;

        % LFP data stream
        if length(openEphys_metadata.continuous)>1
            session.extracellular.srLfp = openEphys_metadata.continuous(2).sample_rate;
        end

        % Events (if present)
        if isfield(openEphys_metadata,'events') && ~isempty(openEphys_metadata.events)
            session.timeSeries.dig.fileName = [openEphys_metadata.events{1}.folder_name,'timestamps.npy'];
            session.timeSeries.dig.fileFormat = 'npy';
            session.timeSeries.dig.precision  = openEphys_metadata.events{1}.type;
            session.timeSeries.dig.nChannels = 1;
            session.timeSeries.dig.sr = openEphys_metadata.events{1}.sample_rate;
            session.timeSeries.dig.equipment = 'Open Ephys Neuropix-PXI';
            session.timeSeries.dig.nSamples = [];
            session.timeSeries.dig.leastSignificantBit = [];
        else
            warning('No events found in the structure.oebin file.');
        end

        % Plot (match your style)
        figure
        plot(chanCoords.x,chanCoords.y,'.k'), hold on
        x = session.extracellular.chanCoords.x;
        y = session.extracellular.chanCoords.y;
        plot(x,y,'s','MarkerFaceColor',[0.9 0.3 0.1])
        xlabel('X position (µm)'), ylabel('Y position (µm)'), title(['Neuropixel 1.0 site selection: ', file1])
        axis equal

    else
        % ================= NP2.0 path (UNCHANGED) =================
        % Electrode groups and channel mapping (keep your current logic)
        channelmapping = zeros(1, session.extracellular.nChannels); % Pre-allocate for efficiency
        for i = 1:session.extracellular.nChannels
            if isstruct(openEphys_metadata.continuous(stream_id).channels) 
                if isfield(openEphys_metadata.continuous(stream_id).channels(i),'channel_metadata')
                    channelmapping(i) = openEphys_metadata.continuous(stream_id).channels(i).channel_metadata.value;
                end
            elseif isfield(openEphys_metadata.continuous(stream_id).channels{i},'channel_metadata')
                channelmapping(i) = openEphys_metadata.continuous(stream_id).channels{i}.channel_metadata.value;
            end
        end

        chanCoords = generateChanCoords_Neuropixel2;

        session.extracellular.chanCoords = chanCoords;
        session.extracellular.chanCoords.x = chanCoords.x(channelmapping + 1); % NP2.0 mapping preserved
        session.extracellular.chanCoords.y = chanCoords.y(channelmapping + 1);
        session.extracellular.chanCoords.shank = chanCoords.shank(channelmapping + 1);
        session.extracellular.chanCoords.channels = chanCoords.channels(channelmapping + 1);

        session.extracellular.electrodeGroups.channels = {};
        
        shanks = unique(session.extracellular.chanCoords.shank);
        for j = 1:length(shanks)
            session.extracellular.electrodeGroups.channels{j} = find(session.extracellular.chanCoords.shank == shanks(j));
            [~,idx] = sortrows([...
                session.extracellular.chanCoords.x(session.extracellular.electrodeGroups.channels{j});...
                session.extracellular.chanCoords.y(session.extracellular.electrodeGroups.channels{j})...
                ]',[2,1],{'descend','ascend'});
            session.extracellular.electrodeGroups.channels{j} = session.extracellular.electrodeGroups.channels{j}(idx);
            session.extracellular.electrodeGroups.label{j} = ['shanks',num2str(shanks(j))]; 
        end
        session.extracellular.nElectrodeGroups = length(shanks);

        figure
        plot(chanCoords.x,chanCoords.y,'.k'), hold on
        site_cmap = hot(session.extracellular.nElectrodeGroups);
        for j = 1:session.extracellular.nElectrodeGroups
            x = session.extracellular.chanCoords.x(session.extracellular.electrodeGroups.channels{j});
            y = session.extracellular.chanCoords.y(session.extracellular.electrodeGroups.channels{j});
            plot(x,y,'s',MarkerFaceColor=site_cmap(j,:))
        end
        xlabel('X position (µm)'), ylabel('Y position (µm)'), title(['Neuropixel site selection: ', file1])
        axis equal

        session.extracellular.spikeGroups.channels = session.extracellular.electrodeGroups.channels;
        session.extracellular.nSpikeGroups = session.extracellular.nElectrodeGroups;

        % LFP data stream
        if length(openEphys_metadata.continuous)>1
            session.extracellular.srLfp = openEphys_metadata.continuous(2).sample_rate;
        end

        % Events (if present)
        if ~isempty(openEphys_metadata.events)
            session.timeSeries.dig.fileName = [openEphys_metadata.events{1}.folder_name,'timestamps.npy'];
            session.timeSeries.dig.fileFormat = 'npy';
            session.timeSeries.dig.precision  = openEphys_metadata.events{1}.type;
            session.timeSeries.dig.nChannels = 1;
            session.timeSeries.dig.sr = openEphys_metadata.events{1}.sample_rate;
            session.timeSeries.dig.equipment = 'Open Ephys Neuropix-PXI';
            session.timeSeries.dig.nSamples = [];
            session.timeSeries.dig.leastSignificantBit = [];
        else
            warning('No events found in the structure.oebin file.');
        end
    end
end
