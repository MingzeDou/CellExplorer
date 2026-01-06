clc;
clear all;

% Load data and process behavior (unchanged)
basepath = 'D:\petersen_lab\data\MV03\MV03_2025-12-21_12-00-27';
[~, basename] = fileparts(basepath);
cd(basepath)

% Load session and tracking data
sessionA = loadSession(basepath,[basename '_ProbeA'],'showGUI',false);
sessionA = loadOpenEphysSettingsFile(fullfile(sessionA.general.basePath, 'structure.oebin'), sessionA, 'probeLetter', 'B');
session = sessionA;
% saveStruct(session);

sessionA = preprocessOpenEphysData('session', session, 'probeLetter', 'B');

% Load digital pulses
openephysDigB = loadOpenEphysDigitalNidaq(sessionA, 'probeLetter', 'B', 'channelNum', 1);

% Load and align behavior session, Probe A, epoch 2
epochToProcess = 2;

% Parameters
scaling_factor = 1;
offset_origin = [5,-5,0];
offset_rigid_body = [5,-5,0];

% Load tracking data
sessionA.behavioralTracking{1}.filenames = 'Take 2025-12-21 12.36.22 PM.csv';
sessionA.behavioralTracking{1}.epoch = epochToProcess;
sessionA.behavioralTracking{1}.framerate = 120;

circular_track_2 = loadOptitrack('session',sessionA,'dataName',['circular_track_' num2str(epochToProcess)],...
    'offset_origin',offset_origin,'scaling_factor',scaling_factor);

% Get on timestamps for the specified epoch
epochOnTimestamps = openephysDigB.on{1,1};

% Truncate TTL timestamps to match behavior frames
epochOnTimestamps = epochOnTimestamps(1:circular_track_2.framesPrFile);

% Report stats
fprintf('\nBehavior Session (Epoch %d):\n', epochToProcess);
fprintf('Behavior frames: %d\n', circular_track_2.framesPrFile);
fprintf('TTL timestamps after truncation: %d\n', length(epochOnTimestamps));
fprintf('Duration of final data: %.2f seconds\n', epochOnTimestamps(end) - epochOnTimestamps(1));

% Update timestamps and metadata
circular_track_2.timestamps = epochOnTimestamps;
circular_track_2.timestamps_reference = 'ephys';
circular_track_2.nSamples = length(epochOnTimestamps);

% Save tracking data
saveStruct(circular_track_2,'behavior','session',sessionA);

% maze parameters
maze = {};
maze.type = 'theta';
maze.radius_in = 96.5/2;
maze.radius_out =  116.5/2;
maze.arm_half_width = 4;
maze.cross_radii = 47.9;
maze.rim_buffer = 10;
maze.polar_rho_limits = [40,75]; % 40,?
maze.polar_theta_limits = [15,2.8*maze.radius_in]; % In units of cm
maze.pos_x_limits = [-13,13]; %
maze.pos_y_limits = [-35,40];

subplot(1,2,1)
if exist('plot_ThetaMaze.m','file')
	plot_ThetaMaze(maze)
end

% Get trials from behavior, Epoch 2
circular_track_2 = getTrials_thetamaze(circular_track_2,maze, 1);

% Linearizing and defining boundaries, Epoch 2
circular_track_2 = linearize_theta_maze(circular_track_2,maze);

% Generating left_right states data
circular_track_2.states.left_right = nan(size(circular_track_2.timestamps));
for i = 1:circular_track_2.trials.alternation.nTrials
    circular_track_2.states.left_right(circular_track_2.trials.alternation.trials==i) = circular_track_2.states.left_right(i);
end
circular_track_2.stateNames.left_right = {'Left','Right'};

% Saving behavioral data
saveStruct(circular_track_2,'behavior','session',session);

% Extract trial parameters
trial_start = circular_track_2.trials.alternation.start;
trial_end = circular_track_2.trials.alternation.end;
nTrials = circular_track_2.trials.alternation.nTrials;

% Determine left vs right trials based on polar_theta values
left_right = zeros(1, nTrials);
for i = 1:nTrials
    % Find indices corresponding to trial timepoints
    idx_start = find(circular_track_2.timestamps >= trial_start(i), 1, 'first');
    idx_end = find(circular_track_2.timestamps <= trial_end(i), 1, 'last');
    
    if isempty(idx_start) || isempty(idx_end) || idx_start >= idx_end
        continue; % Skip problematic trials
    end
    
    % Sample the end of the trial (last 20% of timepoints)
    window_size = ceil(0.2 * (idx_end - idx_start + 1));
    window = (idx_end - window_size + 1) : idx_end;
    
    % Get polar_theta values at the end of the trial
    theta_values = circular_track_2.position.polar_theta(window);
    mean_theta = mean(theta_values, 'omitnan');
    
    % Classify based on polar_theta (following the supervisor's approach)
    % Left = 1, Right = 2 (matches the example plot)
    if mean_theta < -5
        left_right(i) = 1; % Left trial
    elseif mean_theta > 5
        left_right(i) = 2; % Right trial
    else
        % If unclear, try using x coordinate as fallback
        mean_x = mean(circular_track_2.position.x(window), 'omitnan');
        if mean_x < 0
            left_right(i) = 1; % Left trial
        else
            left_right(i) = 2; % Right trial
        end
    end
end

% Store results in trial structure
circular_track_2.trials.alternation.left_right = left_right;

% Set state names
circular_track_2.stateNames.left_right = {'Left', 'Right'};

%% Plot trial times for each trial
figure('Color', 'w', 'Position', [100, 100, 1200, 600]);

% Calculate trial durations
trial_durations = trial_end - trial_start;

% Create subplot for trial times
subplot(2, 1, 1);
bar(1:nTrials, trial_durations, 'FaceColor', [0.4, 0.6, 0.8]);
xlabel('Trial Number');
ylabel('Trial Duration (s)');
title('Trial Duration for Each Trial');
grid on;

% Color code by left/right trials
hold on;
left_trials = find(left_right == 1);
right_trials = find(left_right == 2);
if ~isempty(left_trials)
    bar(left_trials, trial_durations(left_trials), 'FaceColor', [0.8, 0.4, 0.4]);
end
if ~isempty(right_trials)
    bar(right_trials, trial_durations(right_trials), 'FaceColor', [0.4, 0.8, 0.4]);
end
legend({'All Trials', 'Left Trials', 'Right Trials'}, 'Location', 'best');

% Create subplot for trial start times
subplot(2, 1, 2);
trial_start_relative = trial_start - trial_start(1); % Relative to first trial
plot(1:nTrials, trial_start_relative, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.2, 0.4, 0.8]);
xlabel('Trial Number');
ylabel('Trial Start Time (s)');
title('Trial Start Times (Relative to First Trial)');
grid on;

% Color code by left/right trials
hold on;
if ~isempty(left_trials)
    plot(left_trials, trial_start_relative(left_trials), 'o', 'MarkerSize', 8, 'Color', [0.8, 0.4, 0.4]);
end
if ~isempty(right_trials)
    plot(right_trials, trial_start_relative(right_trials), 'o', 'MarkerSize', 8, 'Color', [0.4, 0.8, 0.4]);
end
legend({'All Trials', 'Left Trials', 'Right Trials'}, 'Location', 'best');

%% Plot trajectory for trials
figure('Color', 'w', 'Position', [100, 100, 1200, 800]);

% Plot maze layout
subplot(1, 2, 1);
if exist('plot_ThetaMaze.m','file')
    plot_ThetaMaze(maze);
    hold on;
end

% Plot all trajectories
all_x = circular_track_2.position.x;
all_y = circular_track_2.position.y;
plot(all_x, all_y, 'b-', 'LineWidth', 0.5, 'Alpha', 0.3);
title('All Trajectories');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');
axis equal;
grid on;

% Plot individual trial trajectories
subplot(1, 2, 2);
if exist('plot_ThetaMaze.m','file')
    plot_ThetaMaze(maze);
    hold on;
end

% Use different colors for left and right trials
colors = lines(nTrials);
for i = 1:min(nTrials, 20) % Plot first 20 trials to avoid overcrowding
    % Get trial indices
    idx_start = find(circular_track_2.timestamps >= trial_start(i), 1, 'first');
    idx_end = find(circular_track_2.timestamps <= trial_end(i), 1, 'last');
    
    if isempty(idx_start) || isempty(idx_end) || idx_start >= idx_end
        continue;
    end
    
    % Extract trajectory for this trial
    trial_x = circular_track_2.position.x(idx_start:idx_end);
    trial_y = circular_track_2.position.y(idx_start:idx_end);
    
    % Choose color based on left/right classification
    if left_right(i) == 1 % Left trial
        color = [0.8, 0.4, 0.4]; % Red
    else % Right trial
        color = [0.4, 0.8, 0.4]; % Green
    end
    
    % Plot trajectory with arrows to show direction
    plot(trial_x, trial_y, '-', 'Color', color, 'LineWidth', 2);
    
    % Add arrow at the middle of the trajectory
    mid_idx = round(length(trial_x) / 2);
    if mid_idx > 1 && mid_idx < length(trial_x)
        quiver(trial_x(mid_idx), trial_y(mid_idx), ...
               trial_x(mid_idx+1) - trial_x(mid_idx), ...
               trial_y(mid_idx+1) - trial_y(mid_idx), ...
               0, 'Color', color, 'LineWidth', 2, 'MaxHeadSize', 0.5);
    end
end

title('Individual Trial Trajectories (First 20 Trials)');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');
axis equal;
grid on;

% Add legend for left/right trials
legend({'Maze', 'Left Trials', 'Right Trials'}, 'Location', 'best');

% Create a separate figure for detailed trajectory analysis
figure('Color', 'w', 'Position', [100, 100, 1400, 1000]);

% Plot speed over time for all trials
subplot(3, 1, 1);
plot(circular_track_2.timestamps, circular_track_2.smoothedSpeed, 'b-', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Speed (cm/s)');
title('Speed Over Time');
grid on;

% Mark trial boundaries
for i = 1:nTrials
    line([trial_start(i), trial_start(i)], ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1);
    line([trial_end(i), trial_end(i)], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
end

% Plot position over time
subplot(3, 1, 2);
plot(circular_track_2.timestamps, circular_track_2.position.x, 'r-', 'LineWidth', 1, 'DisplayName', 'X Position');
hold on;
plot(circular_track_2.timestamps, circular_track_2.position.y, 'b-', 'LineWidth', 1, 'DisplayName', 'Y Position');
xlabel('Time (s)');
ylabel('Position (cm)');
title('Position Over Time');
legend('Location', 'best');
grid on;

% Mark trial boundaries
for i = 1:nTrials
    line([trial_start(i), trial_start(i)], ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1);
    line([trial_end(i), trial_end(i)], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
end

% Plot polar coordinates
subplot(3, 1, 3);
plot(circular_track_2.timestamps, circular_track_2.position.polar_theta, 'm-', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Polar Theta (rad)');
title('Polar Angle Over Time');
grid on;

% Mark trial boundaries
for i = 1:nTrials
    line([trial_start(i), trial_start(i)], ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1);
    line([trial_end(i), trial_end(i)], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
end

% Display summary statistics
fprintf('\n=== Trial Summary ===\n');
fprintf('Total number of trials: %d\n', nTrials);
fprintf('Left trials: %d\n', sum(left_right == 1));
fprintf('Right trials: %d\n', sum(left_right == 2));
fprintf('Mean trial duration: %.2f ± %.2f s\n', mean(trial_durations), std(trial_durations));
fprintf('Total recording duration: %.2f s\n', max(circular_track_2.timestamps) - min(circular_track_2.timestamps));