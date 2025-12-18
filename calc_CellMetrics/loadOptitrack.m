function optitrack = loadOptitrack(varargin)
% Loads position tracking data from Optitrack (csv files) to CellExplorer behavior data container
% https://cellexplorer.org/datastructure/data-structure-and-format/#behavior
%
% Example calls
% optitrack = loadOptitrack('session',session)
% optitrack = loadOptitrack('basepath',basepath,'basename',basename,'filenames',filenames)
% theta_maze = loadOptitrack('session',session,'dataName','theta_maze')
% linear_track = loadOptitrack('session',session,'dataName','linear_track')

p = inputParser;
addParameter(p,'session', [], @isstruct); % A session struct
addParameter(p,'basepath', pwd, @isstr); % Basepath of the session
addParameter(p,'basename', [], @isstr); % Name of the session
addParameter(p,'dataName','',@isstr); % Any renaming of the behavior struct
addParameter(p,'filenames', '',@iscellstr); % List of tracking files 
addParameter(p,'scaling_factor', 1, @isnumeric); % A scaling-factor to apply to the x,y,z data.
addParameter(p,'offset_origin', [0,0,0], @isnumeric); % An offset to apply to the x,y,z data (shift of origin). Applied after the scaling. In cm
addParameter(p,'rotation', [], @isnumeric); % A rotation to apply to the x,y data (in degrees). Applied after the offset.
addParameter(p,'saveMat', true, @islogical); % Creates behavior mat file
addParameter(p,'plotFig', true, @islogical); % Creates plot with behavior
addParameter(p,'saveFig', true, @islogical); % Save figure with behavior to summary folder
parse(p,varargin{:})

parameters = p.Results;

if ~isempty(parameters.session)
    session = parameters.session;
    basename = session.general.name;
    basepath = session.general.basePath;
    filenames = session.behavioralTracking{1}.filenames;
else
    basepath = p.Results.basepath;
    basename = p.Results.basename;
end

if isempty(basename)
    basename = basenameFromBasepath(basepath);
end

if isempty(parameters.dataName)
    dataName = 'optitrack';
else
    dataName = parameters.dataName;
end


if ~isempty(parameters.filenames)
    filenames = parameters.filenames;
end

if ischar(filenames)
    filenames = {filenames};
end
% filename = [datapath recording '/' recordings(id).tracking_file];
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';
header_length = 7;


% Loading files
file1 = fullfile(basepath,filenames{1});
if ~exist(file1,'file')
   error('Behavior file does not exist') 
end
fileID = fopen(file1,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false);
fclose(fileID);
FramesPrFile = size(dataArray{1}(header_length:end),1);
if length(filenames) > 1
    for i = 2:length(filenames)
        file0 = fullfile(basepath,filenames{i});
        fileID = fopen(file0,'r');
        dataArray_temp = textscan(fileID, formatSpec, 'Delimiter', ',',  'ReturnOnError', false);
        fclose(fileID);
        for j = 1:length(dataArray)
            dataArray{j} = [dataArray{j};dataArray_temp{j}(header_length:end)];
        end
        FramesPrFile = [FramesPrFile, size(dataArray_temp{1}(header_length:end),1)];
    end
end

optitrack_temp = [];
optitrack_temp.Frame = str2double(dataArray{1}(header_length:end));
optitrack_temp.Time = str2double(dataArray{2}(header_length:end));
optitrack_temp.Xr = str2double(dataArray{3}(header_length:end));
optitrack_temp.Yr = str2double(dataArray{4}(header_length:end));
optitrack_temp.Zr = str2double(dataArray{5}(header_length:end));
optitrack_temp.Wr = str2double(dataArray{6}(header_length:end));
optitrack_temp.X = str2double(dataArray{7}(header_length:end));
optitrack_temp.Y = str2double(dataArray{8}(header_length:end));
optitrack_temp.Z = str2double(dataArray{9}(header_length:end));
optitrack_temp.TotalFrames = str2double(dataArray{12}(1));
optitrack_temp.TotalExportedFrames = str2double(dataArray{14}(1));
optitrack_temp.RotationType = dataArray{16}(1);
optitrack_temp.LenghtUnit = dataArray{18}(1);
optitrack_temp.CoorinateSpace = dataArray{20}(1);

frameRate_old = str2double(dataArray{6}{1});
if isempty(frameRate_old) || isnan(frameRate_old)
    % If old format is invalid, use new format (index 8)
    optitrack_temp.FrameRate = str2double(dataArray{8}{1});
else
    % Use old format value if valid
    optitrack_temp.FrameRate = frameRate_old;
end

clear dataArray
clearvars filename formatSpec fileID dataArray header_length;

%% ========== NaN INTERPOLATION FOR POSITION DATA ==========
fprintf('\n=== Position NaN Detection & Interpolation ===\n');

% Check for NaN values
nan_X = isnan(optitrack_temp.X);
nan_Y = isnan(optitrack_temp.Y);
nan_Z = isnan(optitrack_temp.Z);
nan_any = nan_X | nan_Y | nan_Z;

fprintf('NaN counts before interpolation:\n');
fprintf('  X: %d (%.2f%%)\n', sum(nan_X), 100*sum(nan_X)/length(nan_X));
fprintf('  Y: %d (%.2f%%)\n', sum(nan_Y), 100*sum(nan_Y)/length(nan_Y));
fprintf('  Z: %d (%.2f%%)\n', sum(nan_Z), 100*sum(nan_Z)/length(nan_Z));

if any(nan_any)
    % Create valid data mask
    valid_mask = ~nan_any;
    valid_indices = find(valid_mask);
    
    if sum(valid_mask) > 10  % Need at least 10 valid points to interpolate
        % Interpolate position data using linear interpolation
        optitrack_temp.X = interp1(valid_indices, optitrack_temp.X(valid_mask), ...
            1:length(optitrack_temp.X), 'linear', 'extrap')';
        optitrack_temp.Y = interp1(valid_indices, optitrack_temp.Y(valid_mask), ...
            1:length(optitrack_temp.Y), 'linear', 'extrap')';
        optitrack_temp.Z = interp1(valid_indices, optitrack_temp.Z(valid_mask), ...
            1:length(optitrack_temp.Z), 'linear', 'extrap')';
        
        % For orientation quaternions - interpolate and renormalize
        nan_quat = isnan(optitrack_temp.Xr) | isnan(optitrack_temp.Yr) | ...
                   isnan(optitrack_temp.Zr) | isnan(optitrack_temp.Wr);
        
        if any(nan_quat)
            valid_quat_mask = ~nan_quat;
            valid_quat_indices = find(valid_quat_mask);
            
            if sum(valid_quat_mask) > 10
                optitrack_temp.Xr = interp1(valid_quat_indices, optitrack_temp.Xr(valid_quat_mask), ...
                    1:length(optitrack_temp.Xr), 'linear', 'extrap')';
                optitrack_temp.Yr = interp1(valid_quat_indices, optitrack_temp.Yr(valid_quat_mask), ...
                    1:length(optitrack_temp.Yr), 'linear', 'extrap')';
                optitrack_temp.Zr = interp1(valid_quat_indices, optitrack_temp.Zr(valid_quat_mask), ...
                    1:length(optitrack_temp.Zr), 'linear', 'extrap')';
                optitrack_temp.Wr = interp1(valid_quat_indices, optitrack_temp.Wr(valid_quat_mask), ...
                    1:length(optitrack_temp.Wr), 'linear', 'extrap')';
                
                % Normalize quaternions after interpolation
                quat_norm = sqrt(optitrack_temp.Xr.^2 + optitrack_temp.Yr.^2 + ...
                               optitrack_temp.Zr.^2 + optitrack_temp.Wr.^2);
                optitrack_temp.Xr = optitrack_temp.Xr ./ quat_norm;
                optitrack_temp.Yr = optitrack_temp.Yr ./ quat_norm;
                optitrack_temp.Zr = optitrack_temp.Zr ./ quat_norm;
                optitrack_temp.Wr = optitrack_temp.Wr ./ quat_norm;
                
                fprintf('  ✓ Quaternion orientation interpolated\n');
            end
        end
        
        % Verify interpolation removed NaNs
        remaining_nans = sum(isnan(optitrack_temp.X) | isnan(optitrack_temp.Y) | isnan(optitrack_temp.Z));
        fprintf('✓ Position interpolation complete: %d NaNs remaining\n', remaining_nans);
    else
        warning('Too few valid data points (%d) to interpolate NaNs', sum(valid_mask));
    end
else
    fprintf('✓ No NaNs detected - data is clean\n');
end

%% ========== END NaN INTERPOLATION ==========

% getting position out in cm, and flipping Z and Y axis
position3D = 100*[-optitrack_temp.X,optitrack_temp.Z,optitrack_temp.Y]*parameters.scaling_factor + parameters.offset_origin;

% Rotationg X,Y coordinates around origin if a rotation parameter is given
if ~isempty(parameters.rotation)
    position3D(:,1) = position3D(:,1) * cosd(parameters.rotation) - position3D(:,2) * sind(parameters.rotation);
    position3D(:,2) = position3D(:,1) * sind(parameters.rotation) + position3D(:,2) * cosd(parameters.rotation);
end

% Estimating the speed of the rat
animal_speed = [optitrack_temp.FrameRate*sqrt(sum(diff(position3D)'.^2)),0];
animal_speed = nanconv(animal_speed,ones(1,10)/10,'edge');  % Original smoothing

% Additional smoothing with larger window
window_size = 25;
smoothed_speed = nanconv(animal_speed, ones(1,window_size)/window_size, 'edge');
animal_acceleration = [0,optitrack_temp.FrameRate*diff(smoothed_speed)];

% Adding output struct
optitrack_temp.position3D = position3D';

% Generating buzcode fields and output struct
optitrack.timestamps = optitrack_temp.Time(:)';
optitrack.timestamps_reference = 'optitrack';
optitrack.sr = optitrack_temp.FrameRate;
optitrack.position.x = optitrack_temp.position3D(1,:);
optitrack.position.y = optitrack_temp.position3D(2,:);
optitrack.position.z = optitrack_temp.position3D(3,:);
optitrack.position.units = 'centimeters';
optitrack.position.referenceFrame = 'global';
optitrack.position.coordinateSystem = 'cartesian';
optitrack.speed = animal_speed;  % Original speed
optitrack.smoothedSpeed = smoothed_speed;  % Added smoothed speed
optitrack.acceleration = animal_acceleration;
optitrack.orientation.x = -optitrack_temp.Xr;
optitrack.orientation.y = optitrack_temp.Zr;
optitrack.orientation.z = optitrack_temp.Yr;
optitrack.orientation.w = optitrack_temp.Wr;
% Create rotation matrix from quaternion components
R11 = 1 - 2*(optitrack.orientation.y.^2 + optitrack.orientation.z.^2);
R12 = 2*(optitrack.orientation.x.*optitrack.orientation.y - optitrack.orientation.w.*optitrack.orientation.z);
R21 = 2*(optitrack.orientation.x.*optitrack.orientation.y + optitrack.orientation.w.*optitrack.orientation.z);
R22 = 1 - 2*(optitrack.orientation.x.^2 + optitrack.orientation.z.^2);

% Calculate forward direction vector (use first column of rotation matrix)
% This represents where the animal is facing
forward_x = R11;
forward_y = R21;

% Calculate heading as angle in the XY plane
heading = atan2(forward_y, forward_x);

% Rotate heading 90 degrees counterclockwise (add pi/2)
heading_rotated = heading + (pi/2);

% Wrap to pi range
optitrack.orientation.heading = wrapToPi(heading_rotated);

optitrack.orientation.rotationType = optitrack_temp.RotationType;
optitrack.nSamples = numel(optitrack.timestamps);

if exist('FramesPrFile')
    optitrack.framesPrFile = FramesPrFile;
end

% Attaching info about how the data was processed
optitrack.processinginfo.function = 'loadOptitrack';
optitrack.processinginfo.version = 1;
optitrack.processinginfo.date = now;
optitrack.processinginfo.params.basepath = basepath;
optitrack.processinginfo.params.basename = basename;
optitrack.processinginfo.params.filenames = filenames;
optitrack.processinginfo.params.dataName = dataName;
optitrack.processinginfo.params.rotation = parameters.rotation;
optitrack.processinginfo.params.scaling_factor = parameters.scaling_factor;
optitrack.processinginfo.params.offset_origin = parameters.offset_origin;

try
    optitrack.processinginfo.username = char(java.lang.System.getProperty('user.name'));
    optitrack.processinginfo.hostname = char(java.net.InetAddress.getLocalHost.getHostName);
catch
    disp('Failed to retrieve system info.')
end

% Saving data
if parameters.saveMat
    saveStruct(optitrack,'behavior','dataName',dataName,'session',session);
end

% Plotting
if parameters.plotFig
    fig1 = figure;
    subplot(1,2,1)
    % plot3(position3D(:,1),position3D(:,2),position3D(:,3)), title('Position'), xlabel('X (cm)'), ylabel('Y (cm)'), zlabel('Z (cm)'),axis tight, view(2), hold on
    plot3(position3D(:,1),position3D(:,2),position3D(:,3)), title('Position'), xlabel('X (cm)'), ylabel('Y (cm)'), zlabel('Z (cm)'),axis tight, hold on
    subplot(1,2,2)
    plot3(position3D(:,1),position3D(:,2),smoothed_speed), hold on
    xlabel('X (cm)'), ylabel('Y (cm)'),zlabel('Speed (cm/s)'), axis tight
    
    % Saving a summary figure for all cells
    if parameters.saveFig
        timestamp = datestr(now, '_dd-mm-yyyy_HH.MM.SS');
        ce_savefigure(fig1,basepath,[basename, '.',dataName,'.behavior' timestamp])
        disp(['loadOptitrack: Summary figure saved to ', fullfile(basepath, 'SummaryFigures', [basename, '.', dataName, '.behavior', timestamp]),'.png'])
    end
end
