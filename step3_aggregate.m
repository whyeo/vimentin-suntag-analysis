% step3_aggregate.m
% -------------------------------------------------------------------------
% Title:    step3_aggregate
% Purpose:  Aggregate per-file analysis results produced by step2, pad and
%           concatenate track-level matrices across files, and produce
%           condition-level summary figures for visualization.
%
% Usage:    Edit the 'experiment_date' and 'experiment_conditions' variables
%           below and run this script from the project root.
%
% Inputs:   graphs/<experiment_date>/<condition>/*_MSD_DP.mat
% Outputs:  graphs/<experiment_date>/<condition>/track_analysis_*.png/.eps
%
% Notes:    - Verify 'pixel_size_nm', 'frame_rate_hz' and normalization constants
%            to ensure correct physical units for distances and speeds.
%
% Author:   Wei Hong Yeo
% Date:     2026-01-04
% -------------------------------------------------------------------------

%% User parameters (replace values as needed)
% Set the experiment date (folder name under examples/)
experiment_date = '2023-12-17';
% list of experimental condition folders under examples/<date>/
experiment_conditions = {'Condition1', 'Condition2'};

% Imaging and analysis parameters (verify for your setup)
pixel_size_nm = 130;    % camera pixel size in nanometers
photon_normalization = 500; % intensity normalization factor (photons)
frame_rate_hz = 2;      % imaging frame rate in Hz
frame_interval_s = 1 / frame_rate_hz;
plot_tracks = true;     % toggle drawing individual tracks in diagnostics

% optional color helper
try
    colors_accent = colors('accent');
catch
    colors_accent = [0.9,0.3,0.2; 0.2,0.6,0.2; 0.3,0.5,0.9]; % fallback
end

% track type codes
STATIONARY_TRACK = 0;
SHORT_TRACK = 1;
LONG_TRACK = 2;

% ensure top-level output folder
ensureDir(fullfile('graphs', experiment_date));

n_conditions = numel(experiment_conditions);
for i_condition = 1:n_conditions
    cond_name = experiment_conditions{i_condition};
    ts_folder = fullfile('examples', experiment_date, cond_name, 'tsoutput');
    filelist = dir(fullfile(ts_folder, '*-sub.csv'));
    n_files = numel(filelist);

    % initialize aggregated containers
    aggregated_track_mean_squared_displacement = [];
    aggregated_track_alpha = [];
    aggregated_track_directional_ratios = [];
    aggregated_track_directional_autocorrelation = [];
    aggregated_track_fractional_long_distance = [];
    aggregated_track_fractional_short_distance = [];
    aggregated_track_fractional_stationary = [];
    aggregated_track_type = [];
    aggregated_track_mean_speed = [];
    aggregated_track_maximum_velocity = [];

    % ensure condition output folder exists
    cond_out = fullfile('graphs', experiment_date, cond_name);
    ensureDir(cond_out);

    for i_file = 1:n_files
        csv_name = filelist(i_file).name;
        save_data_file = fullfile(cond_out, sprintf('%s_MSD_DP.mat', csv_name(1:end-8)));
        if ~exist(save_data_file, 'file')
            warning('  Missing per-file MAT: %s (skipping)', save_data_file);
            continue;
        end

        S = load(save_data_file);
        % expected variables inside S: track_mean_squared_displacement, track_alpha,
        % track_directional_ratios, track_directional_autocorrelation, timelags,
        % track_type, track_fractional_long_distance, track_fractional_short_distance,
        % track_fractional_stationary, track_mean_speed, track_maximum_velocity

        % make sure matrices have compatible number of columns (lags). Pad as needed.
        % pad MSD
        if isfield(S, 'track_mean_squared_displacement')
            [aggregated_track_mean_squared_displacement, S.track_mean_squared_displacement] = padToSameCols(aggregated_track_mean_squared_displacement, S.track_mean_squared_displacement);
            aggregated_track_mean_squared_displacement = [aggregated_track_mean_squared_displacement; S.track_mean_squared_displacement];
        end

        % pad directional ratios and autocorr
        if isfield(S, 'track_directional_ratios')
            [aggregated_track_directional_ratios, S.track_directional_ratios] = padToSameCols(aggregated_track_directional_ratios, S.track_directional_ratios);
            aggregated_track_directional_ratios = [aggregated_track_directional_ratios; S.track_directional_ratios];
        end
        if isfield(S, 'track_directional_autocorrelation')
            [aggregated_track_directional_autocorrelation, S.track_directional_autocorrelation] = padToSameCols(aggregated_track_directional_autocorrelation, S.track_directional_autocorrelation);
            aggregated_track_directional_autocorrelation = [aggregated_track_directional_autocorrelation; S.track_directional_autocorrelation];
        end  
        aggregated_track_mean_squared_displacement = [aggregated_track_mean_squared_displacement; S.track_mean_squared_displacement];
        aggregated_track_alpha = [aggregated_track_alpha; S.track_alpha];
        aggregated_track_directional_ratios = [aggregated_track_directional_ratios; S.track_directional_ratios];
        aggregated_track_directional_autocorrelation = [aggregated_track_directional_autocorrelation; S.track_directional_autocorrelation];
        aggregated_track_fractional_long_distance = [aggregated_track_fractional_long_distance; S.track_fractional_long_distance];
        aggregated_track_fractional_short_distance = [aggregated_track_fractional_short_distance; S.track_fractional_short_distance];
        aggregated_track_fractional_stationary = [aggregated_track_fractional_stationary; S.track_fractional_stationary];
        aggregated_track_type = [aggregated_track_type; S.track_type];
        aggregated_track_mean_speed = [aggregated_track_mean_speed; S.track_mean_speed];
        aggregated_track_maximum_velocity = [aggregated_track_maximum_velocity; S.track_maximum_velocity];
    end % per-file loop

    % If no data aggregated for this condition, skip plotting
    if isempty(aggregated_track_type)
        warning('No aggregated data for condition %s - skipping plots.', cond_name);
        continue;
    end

    % construct timelags based on aggregated MSD columns
    n_lags = size(aggregated_track_mean_squared_displacement, 2);
    timelags = (0:n_lags-1) * frame_interval_s;

    % create summary figure (similar layout to original)
    f3 = figure('Name', sprintf('Track Analysis - %s', cond_name), 'NumberTitle','off'); clf;
    f3.Position(3:4) = [900, 600];
    tiledlayout(4,2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % MSD per-class
    nexttile; hold on; box on;
    xlabel('Time lag (s)'); ylabel('Mean Squared Displacement (\\mum^2)');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    try
        y_long = mean(aggregated_track_mean_squared_displacement(aggregated_track_type == LONG_TRACK, :), 1, 'omitnan');
        e_long = std(aggregated_track_mean_squared_displacement(aggregated_track_type == LONG_TRACK, :), 0, 1, 'omitnan');
        errorbar(timelags, y_long, e_long, 'Color', colors_accent(3,:));
    catch
    end
    try
        y_short = mean(aggregated_track_mean_squared_displacement(aggregated_track_type == SHORT_TRACK, :), 1, 'omitnan');
        e_short = std(aggregated_track_mean_squared_displacement(aggregated_track_type == SHORT_TRACK, :), 0, 1, 'omitnan');
        errorbar(timelags, y_short, e_short, 'Color', colors_accent(1,:));
    catch
    end
    try
        y_stat = mean(aggregated_track_mean_squared_displacement(aggregated_track_type == STATIONARY_TRACK, :), 1, 'omitnan');
        e_stat = std(aggregated_track_mean_squared_displacement(aggregated_track_type == STATIONARY_TRACK, :), 0, 1, 'omitnan');
        errorbar(timelags, y_stat, e_stat, 'Color', [0.5,0.5,0.5]);
    catch
    end

    % Directional persistence
    nexttile; hold on; box on; xlabel('Time lag (s)'); ylabel('Directional Persistence (d/D)');
    try
        errorbar(timelags, mean(aggregated_track_directional_ratios(aggregated_track_type == LONG_TRACK, :), 1, 'omitnan'), std(aggregated_track_directional_ratios(aggregated_track_type == LONG_TRACK, :), 0, 1, 'omitnan'), 'Color', colors_accent(3,:));
        errorbar(timelags, mean(aggregated_track_directional_ratios(aggregated_track_type == SHORT_TRACK, :), 1, 'omitnan'), std(aggregated_track_directional_ratios(aggregated_track_type == SHORT_TRACK, :), 0, 1, 'omitnan'), 'Color', colors_accent(1,:));
        errorbar(timelags, mean(aggregated_track_directional_ratios(aggregated_track_type == STATIONARY_TRACK, :), 1, 'omitnan'), std(aggregated_track_directional_ratios(aggregated_track_type == STATIONARY_TRACK, :), 0, 1, 'omitnan'), 'Color', [0.5,0.5,0.5]);
    catch
    end

    % Mean speed (classified)
    nexttile; hold on; box on; xlabel('Mean Speed (\mum/s)'); ylabel('Probability');
    bins = 0:0.02:0.6;
    try
        histogram(aggregated_track_mean_speed(aggregated_track_type == LONG_TRACK), bins, 'Normalization','probability', 'DisplayStyle','stairs', 'EdgeColor', colors_accent(3,:));
        histogram(aggregated_track_mean_speed(aggregated_track_type == SHORT_TRACK), bins, 'Normalization','probability', 'DisplayStyle','stairs', 'EdgeColor', colors_accent(1,:));
        histogram(aggregated_track_mean_speed(aggregated_track_type == STATIONARY_TRACK), bins, 'Normalization','probability', 'DisplayStyle','stairs', 'EdgeColor', [0.5,0.5,0.5]);
    catch
    end

    % Maximum velocity (classified)
    nexttile; hold on; box on; xlabel('Maximum Velocity (\mum/s)'); ylabel('Probability');
    bins2 = 0:0.04:1.2;
    try
        histogram(aggregated_track_maximum_velocity(aggregated_track_type == LONG_TRACK), bins2, 'Normalization','probability', 'DisplayStyle','stairs', 'EdgeColor', colors_accent(3,:));
        histogram(aggregated_track_maximum_velocity(aggregated_track_type == SHORT_TRACK), bins2, 'Normalization','probability', 'DisplayStyle','stairs', 'EdgeColor', colors_accent(1,:));
        histogram(aggregated_track_maximum_velocity(aggregated_track_type == STATIONARY_TRACK), bins2, 'Normalization','probability', 'DisplayStyle','stairs', 'EdgeColor', [0.5,0.5,0.5]);
    catch
    end

    % overall mean speed
    nexttile; hold on; box on; xlabel('Mean Speed (\mum/s)'); ylabel('Probability');
    try
        histogram(aggregated_track_mean_speed, bins, 'Normalization','probability', 'DisplayStyle','stairs', 'EdgeColor', 'k');
    catch
    end

    % overall maximum velocity
    nexttile; hold on; box on; xlabel('Maximum Velocity (\mum/s)'); ylabel('Probability');
    try
        histogram(aggregated_track_maximum_velocity, bins2, 'Normalization','probability', 'DisplayStyle','stairs', 'EdgeColor', 'k');
    catch
    end

    % alpha distribution
    nexttile; hold on; box on; xlabel('Alpha'); ylabel('Probability');
    try
        histogram(aggregated_track_alpha, 0:0.02:1, 'Normalization','probability', 'DisplayStyle','stairs', 'EdgeColor', 'k');
    catch
    end

    % title and save
    sgtitle(sprintf('Track Analysis — %s', cond_name), 'Interpreter', 'none');
    figure_savefile = fullfile(cond_out, ['track_analysis_' cond_name]);
    try
        formatfig(f3, 'publication', [figure_savefile, '.png']);
        formatfig(f3, 'publication', [figure_savefile, '.eps'], 'painters');
    catch
        saveas(f3, [figure_savefile, '.png']);
    end
    close(f3);
end

%% Local helper: padToSameCols must be at end of the script (MATLAB requirement)
function [A_padded, B_padded] = padToSameCols(A, B)
    colsA = size(A,2);
    colsB = size(B,2);
    nmax = max(colsA, colsB);
    if isempty(A)
        A_padded = nan(0, nmax);
    else
        A_padded = A;
        if colsA < nmax
            A_padded(:, colsA+1:nmax) = nan(size(A_padded,1), nmax-colsA);
        end
    end
    if isempty(B)
        B_padded = nan(0, nmax);
    else
        B_padded = B;
        if colsB < nmax
            B_padded(:, colsB+1:nmax) = nan(size(B_padded,1), nmax-colsB);
        end
    end
end

%% Local helper
function ensureDir(folder)
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end