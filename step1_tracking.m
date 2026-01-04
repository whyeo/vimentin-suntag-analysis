% step1_tracking.m
% -------------------------------------------------------------------------
% Title:    step1_tracking
% Purpose:  Load ThunderSTORM CSV localizations, perform tracking and compute
%           per-track statistics (length, displacement, total distance,
%           velocity, MSD). Save per-file caches and diagnostic outputs.
%
% Usage:    Edit the 'experiment_date' and 'experiment_conditions' variables
%           below and run this script from the project root.
%
% Inputs:   examples/<experiment_date>/<condition>/tsoutput/*-sub.csv
% Outputs:  graphs/<experiment_date>/<condition>/* (PNG, XLSX, .mat caches)
%
% Notes:    - Verify 'pixel_size_nm', 'frame_rate_hz' and normalization constants
%            to ensure correct physical units for distances and speeds.
%           - Requires a 'track.m' implementation on the MATLAB path.
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

%% Preallocate containers for assembled results
n_conditions = length(experiment_conditions);
assemble_track_full = cell(1,n_conditions);
assemble_track_lengths = cell(1,n_conditions);
assemble_track_velocities = cell(1,n_conditions);
assemble_track_displacements = cell(1,n_conditions);
assemble_track_mean_speed = cell(1,n_conditions);
assemble_track_total_distance = cell(1,n_conditions);
assemble_track_maximum_distance = cell(1,n_conditions);
assemble_track_mean_squared_displacement = cell(1,n_conditions);
assemble_track_fractional_long_distance = cell(1,n_conditions);
assemble_track_summary_counts = cell(1,n_conditions);

% Ensure top-level output directory exists
out_base = fullfile('graphs', experiment_date);
ensureDir(out_base);

%% Process each experimental condition
for ic = 1:n_conditions
    cond_name = experiment_conditions{ic}; % use {} to extract char from cell
    fprintf('Processing condition: %s\n', cond_name);

    ts_dir = fullfile('examples', experiment_date, cond_name, 'tsoutput');
    file_patterns = fullfile(ts_dir, '*-sub.csv');
    files = dir(file_patterns);

    n_files = length(files);

    % Prepare per-condition containers
    assemble_track_full{ic} = cell(1, n_files);
    assemble_track_lengths{ic} = cell(1, n_files);
    assemble_track_velocities{ic} = cell(1, n_files);
    assemble_track_displacements{ic} = cell(1, n_files);
    assemble_track_mean_speed{ic} = cell(1, n_files);
    assemble_track_total_distance{ic} = cell(1, n_files);
    assemble_track_maximum_distance{ic} = cell(1, n_files);
    assemble_track_mean_squared_displacement{ic} = cell(1, n_files);
    assemble_track_fractional_long_distance{ic} = zeros(1, n_files);
    assemble_track_summary_counts{ic} = zeros(1, n_files);

    % Create folder for this condition's graphs
    cond_out = fullfile(out_base, cond_name);
    ensureDir(cond_out);

    for ifile = 1:n_files
        csv_name = files(ifile).name;
        csv_path = fullfile(ts_dir, csv_name);
        fprintf('  File (%d/%d): %s\n', ifile, n_files, csv_name);

        try
            % Read the ThunderSTORM CSV
            data = readtable(csv_path, 'PreserveVariableNames', true);

            % Extract variables and convert units where appropriate
            x_px = data{:,'x [nm]'} / pixel_size_nm; % convert nm->pixels
            y_px = data{:,'y [nm]'} / pixel_size_nm;
            frames = data{:,'frame'}; % frame indices (integers)
            intensity_norm = data{:,'intensity [photon]'} / photon_normalization;

            % Prepare data for the tracker: [x y intensity frame]
            data_for_tracking = [x_px, y_px, intensity_norm, frames];

            % Tracking parameters for track.m
            tracking_params.mem = 2;
            tracking_params.dim = 3;
            tracking_params.good = 5;
            tracking_params.quiet = 1;

            % Run tracking (track.m must be available on path)
            tracks = track(data_for_tracking, 4, tracking_params);

            % Save track results for this file (optional cache)
            track_cache_name = fullfile(ts_dir, sprintf('%s_track.mat', csv_name(1:end-8)));
            try
                save(track_cache_name, 'tracks');
            catch
                warning('  Could not save track cache: %s', track_cache_name);
            end

            % Collect per-track statistics
            n_tracks = max(tracks(:,end));
            track_full = cell(n_tracks, 1);
            track_length_s = zeros(n_tracks,1);
            track_velocities_um_s = zeros(n_tracks,1);
            track_mean_speed_um_s = zeros(n_tracks,1);
            track_displacement_um = zeros(n_tracks,1);
            track_total_distance_um = zeros(n_tracks,1);
            track_maximum_distance_um = zeros(n_tracks,1);
            track_msd = zeros(n_tracks,1);

            % Create figure for diagnostics
            f = figure(1); clf;
            f.Position(3:4) = [1800,900];

            for it = 1:n_tracks
                idx = tracks(:,end) == it;
                curr = tracks(idx, :); % columns: x, y, intensity, frame, ... track id
                track_full{it} = curr;

                % Length in frames, displacement and distances in pixels
                length_frames = max(curr(:,end-1)) - min(curr(:,end-1));
                displacement_px = norm(curr(end,1:2) - curr(1,1:2));
                total_dist_px = sum(sqrt(sum(diff(curr(:,1:2)).^2,2)));

                % velocity and mean speed in pixels/frame
                velocity_px_per_frame = displacement_px / length_frames;
                mean_speed_px_per_frame = total_dist_px / length_frames;

                % maximum single-step distance (in pixels)
                max_step_px = max(sqrt(sum(diff(curr(:,1:2)).^2,2)));

                % mean squared displacement (MSD) estimate (pixels^2)
                msd_px2 = mean((sqrt(sum(diff(curr(:,1:2)).^2,2))).^2);

                % Convert to physical units (um and um/s)
                length_s = length_frames * frame_interval_s;
                velocity_um_s = velocity_px_per_frame * (pixel_size_nm/1000) / frame_interval_s;
                mean_speed_um_s = mean_speed_px_per_frame * (pixel_size_nm/1000) / frame_interval_s;
                displacement_um = displacement_px * (pixel_size_nm/1000);
                total_dist_um = total_dist_px * (pixel_size_nm/1000);
                max_step_um = max_step_px * (pixel_size_nm/1000);
                msd_um2 = msd_px2 * (pixel_size_nm/1000)^2;

                % Store results
                track_length_s(it) = length_s;
                track_velocities_um_s(it) = velocity_um_s;
                track_mean_speed_um_s(it) = mean_speed_um_s;
                track_displacement_um(it) = displacement_um;
                track_total_distance_um(it) = total_dist_um;
                track_maximum_distance_um(it) = max_step_um;
                track_msd(it) = msd_um2;

                % Optionally plot tracks on subplots
                if plot_tracks
                    subplot(3,3,[1,4]); hold on;
                    if track_maximum_distance_um(it) > 1
                        plot(curr(:,1), curr(:,2), '-', 'LineWidth', 2, 'Color', '#5C78DE');
                        plot(curr(1,1), curr(1,2), 'o', 'MarkerSize', 4, 'MarkerFaceColor', '#FFEC3F', 'MarkerEdgeColor', '#FFEC3F');
                    else
                        plot(curr(:,1), curr(:,2), '-', 'LineWidth', 2, 'Color', '#E6194B');
                    end

                    subplot(3,3,[2,5]); hold on;
                    if track_total_distance_um(it) > 3
                        plot(curr(:,1), curr(:,2), '-', 'LineWidth', 2, 'Color', '#5C78DE');
                        plot(curr(1,1), curr(1,2), 'o', 'MarkerSize', 4, 'MarkerFaceColor', '#FFEC3F', 'MarkerEdgeColor', '#FFEC3F');
                    else
                        plot(curr(:,1), curr(:,2), '-', 'LineWidth', 2, 'Color', '#E6194B');
                    end

                    subplot(3,3,[3,6]); hold on;
                    if track_mean_speed_um_s(it) > 0.1
                        plot(curr(:,1), curr(:,2), '-', 'LineWidth', 2, 'Color', '#5C78DE');
                        plot(curr(1,1), curr(1,2), 'o', 'MarkerSize', 4, 'MarkerFaceColor', '#FFEC3F', 'MarkerEdgeColor', '#FFEC3F');
                    else
                        plot(curr(:,1), curr(:,2), '-', 'LineWidth', 2, 'Color', '#E6194B');
                    end
                end
            end % per-track loop

            % Histograms for the file (maximum distance, total distance, mean speed)
            subplot(3,3,7); hold on;
            histogram(track_maximum_distance_um(track_maximum_distance_um <= 1), 100, 'BinLimits', [0,2], 'BinWidth', 0.1, 'FaceColor', '#E6194B');
            histogram(track_maximum_distance_um(track_maximum_distance_um > 1), 100, 'BinLimits', [0,2], 'BinWidth', 0.1, 'FaceColor', '#5C78DE');
            xlabel('maximum distance (um)'); ylabel('count'); xline(1, '--', 'LineWidth', 2, 'Color', '#000000');

            subplot(3,3,8); hold on;
            histogram(track_total_distance_um(track_total_distance_um <= 3), 100, 'BinLimits', [0,5], 'BinWidth', 0.1, 'FaceColor', '#E6194B');
            histogram(track_total_distance_um(track_total_distance_um > 3), 100, 'BinLimits', [0,5], 'BinWidth', 0.1, 'FaceColor', '#5C78DE');
            xlabel('total distance (um)'); ylabel('count'); xline(3, '--', 'LineWidth',2,'Color','#000000');

            subplot(3,3,9); hold on;
            histogram(track_mean_speed_um_s(track_mean_speed_um_s <= 0.1), 100, 'BinLimits', [0,0.2], 'BinWidth', 0.005, 'FaceColor', '#E6194B');
            histogram(track_mean_speed_um_s(track_mean_speed_um_s > 0.1), 100, 'BinLimits', [0,0.2], 'BinWidth', 0.005, 'FaceColor', '#5C78DE');
            xlabel('mean speed (um/s)'); ylabel('count'); xline(0.1, '--', 'LineWidth',2,'Color','#000000');

            % Fraction of tracks longer than 1 um
            fractional_long_distance = sum(track_maximum_distance_um > 1) / length(track_maximum_distance_um);
            total_tracks_in_file = length(track_maximum_distance_um);

            % If plotting, adjust axes and save figure
            if plot_tracks
                subplot(3,3,[1,4]); axis equal; set(gca, 'YDir', 'reverse');
                subplot(3,3,[2,5]); axis equal; set(gca, 'YDir', 'reverse');
                subplot(3,3,[3,6]); axis equal; set(gca, 'YDir', 'reverse');

                out_png = fullfile(cond_out, sprintf('%s.png', csv_name(1:end-4)));
                try
                    formatfig(f, 'presentation', out_png);
                catch
                    % If formatfig not available, fall back to saveas
                    saveas(f, out_png);
                end
            end

            % Store assembled per-file results
            assemble_track_full{ic}{ifile} = track_full;
            assemble_track_lengths{ic}{ifile} = track_length_s;
            assemble_track_velocities{ic}{ifile} = track_velocities_um_s;
            assemble_track_mean_speed{ic}{ifile} = track_mean_speed_um_s;
            assemble_track_displacements{ic}{ifile} = track_displacement_um;
            assemble_track_total_distance{ic}{ifile} = track_total_distance_um;
            assemble_track_maximum_distance{ic}{ifile} = track_maximum_distance_um;
            assemble_track_mean_squared_displacement{ic}{ifile} = track_msd;

            assemble_track_fractional_long_distance{ic}(ifile) = fractional_long_distance;
            assemble_track_summary_counts{ic}(ifile) = total_tracks_in_file;

            % Save track data and summary to an Excel file (two sheets)
            out_xlsx = fullfile(cond_out, sprintf('%s.xlsx', csv_name(1:end-4)));
            track_data_table = table(tracks(:,1), tracks(:,2), tracks(:,end-1), tracks(:,end), 'VariableNames', {'x','y','frame','track_id'});
            track_statistics_table = table((1:n_tracks)', track_length_s, track_displacement_um, track_total_distance_um, track_velocities_um_s, track_maximum_distance_um, track_msd, 'VariableNames', {'track_id','length_s','displacement_um','total_distance_um','velocity_um_s','maximum_distance_um','msd_um2'});
            try
                writetable(track_data_table, out_xlsx, 'Sheet', 'track_data');
                writetable(track_statistics_table, out_xlsx, 'Sheet', 'track_statistics');
            catch ME
                warning('  Could not write XLSX for %s: %s', csv_name, ME.message);
            end

            fprintf('  Finished processing %s (tracks: %d)\n', csv_name(1:end-4), total_tracks_in_file);

        catch ME
            warning('  Error processing %s: %s', csv_name, ME.message);
            continue; % continue with next file
        end
    end % per-file loop
end % per-condition loop

%% Assemble and plot summary across all files and conditions (example figures)
% Combine per-file results into per-condition vectors/matrices and create
% simple summary figures suitable for sharing. This section is a cleaned and
% minimal refactor of the original plotting/assembly logic.

n_conditions = length(experiment_conditions);
assemble_colors = ["#E6194B","#3CB44B","#5C78DE","#FA9F16","#A366EE"];

% Quick safety: if no data were processed, skip the assembly
if all(cellfun(@isempty, assemble_track_lengths))
    warning('No track data found in assemble arrays — skipping summary assembly.');
else
    % Figure: per-file histograms (maximum distance, mean speed, total distance)
    f6 = figure('Name','Per-file histograms','NumberTitle','off'); clf;
    f6.Position(3:4) = [1200,800];

    % Count number of files total to layout subplots
    total_files = sum(cellfun(@length, assemble_track_lengths));
    if total_files == 0
        warning('No files to plot in per-file histograms.');
    else
        i_subplot = 1;
        for i_condition = 1:n_conditions
            for i_file = 1:length(assemble_track_maximum_distance{i_condition})
                % extract vectors for this file (already in physical units: um, um/s)
                max_dist = assemble_track_maximum_distance{i_condition}{i_file};
                mean_speed = assemble_track_mean_speed{i_condition}{i_file};
                total_dist = assemble_track_total_distance{i_condition}{i_file};

                subplot(total_files,3,3*i_subplot-2); hold on;
                histogram(max_dist, 50, 'BinLimits', [0, max(1.5, prctile(max_dist,95))], 'Normalization','probability', 'FaceColor', assemble_colors(min(i_condition,end)));
                xlim([0, max(1.5, prctile(max_dist,95))]);
                title(sprintf('%s - file %d - max dist (um)', experiment_conditions{i_condition}, i_file));

                subplot(total_files,3,3*i_subplot-1); hold on;
                histogram(mean_speed, 50, 'BinLimits', [0, max(0.2, prctile(mean_speed,95))], 'Normalization','probability', 'FaceColor', assemble_colors(min(i_condition,end)));
                xlim([0, max(0.2, prctile(mean_speed,95))]);
                title('mean speed (um/s)');

                subplot(total_files,3,3*i_subplot); hold on;
                histogram(total_dist, 50, 'BinLimits', [0, max(5, prctile(total_dist,95))], 'Normalization','probability', 'FaceColor', assemble_colors(min(i_condition,end)));
                xlim([0, max(5, prctile(total_dist,95))]);
                title('total distance (um)');

                i_subplot = i_subplot + 1;
            end
        end

        % Save figure (use formatfig if available)
        out_png6 = fullfile('graphs', experiment_date, 'track_maximum_distance.png');
        try
            formatfig(f6, 'presentation', out_png6);
        catch
            saveas(f6, out_png6);
        end
    end

    % Figure: combined per-condition summary plots (flattened across files)
    f8 = figure('Name','Condition summaries','NumberTitle','off'); clf; hold on;
    f8.Position(3:4) = [1600,800];

    % For each condition, concatenate the per-file vectors
    for i_condition = 1:n_conditions
        % flatten vectors (if some files are empty, skip them)
        concat_lengths = vertcat(assemble_track_lengths{i_condition}{:});
        concat_velocities = vertcat(assemble_track_velocities{i_condition}{:});
        concat_mean_speed = vertcat(assemble_track_mean_speed{i_condition}{:});
        concat_displacements = vertcat(assemble_track_displacements{i_condition}{:});
        concat_total_distance = vertcat(assemble_track_total_distance{i_condition}{:});
        concat_maximum_distance = vertcat(assemble_track_maximum_distance{i_condition}{:});
        concat_msd = vertcat(assemble_track_mean_squared_displacement{i_condition}{:});

        % Some conditions may be empty — skip plotting those
        if isempty(concat_lengths)
            continue;
        end

        n_subplot_columns = 7;
        % row index for this condition
        row_idx = i_condition;

        % Helper to compute subplot index in a grid of n_conditions x n_subplot_columns
        base = (row_idx-1)*n_subplot_columns;

        subplot(n_conditions,n_subplot_columns,base+1); hold on;
        histogram(concat_lengths, 25, 'BinLimits', [0,60], 'Normalization','probability', 'FaceColor', assemble_colors(min(i_condition,end)));
        xlim([0,60]); xlabel('track length (s)'); ylabel('probability');

        subplot(n_conditions,n_subplot_columns,base+2); hold on;
        histogram(concat_displacements, 25, 'BinLimits', [0,1], 'Normalization','probability', 'FaceColor', assemble_colors(min(i_condition,end)));
        xlim([0,1]); xlabel('track displacement (um)');

        subplot(n_conditions,n_subplot_columns,base+3); hold on;
        histogram(concat_msd, 25, 'BinLimits', [0,0.05], 'Normalization','probability', 'FaceColor', assemble_colors(min(i_condition,end)));
        xlim([0,0.05]); xlabel('mean squared displacement (um^2)');

        subplot(n_conditions,n_subplot_columns,base+4); hold on;
        histogram(concat_velocities, 25, 'BinLimits', [0,0.2], 'Normalization','probability', 'FaceColor', assemble_colors(min(i_condition,end)));
        xlim([0,0.2]); xlabel('track velocity (um/s)');

        subplot(n_conditions,n_subplot_columns,base+5); hold on;
        histogram(concat_maximum_distance, 25, 'BinLimits', [0,2], 'Normalization','probability', 'FaceColor', assemble_colors(min(i_condition,end)));
        xlim([0,2]); xlabel('maximum distance (um)');

        subplot(n_conditions,n_subplot_columns,base+6); hold on;
        histogram(concat_mean_speed, 25, 'BinLimits', [0,0.2], 'Normalization','probability', 'FaceColor', assemble_colors(min(i_condition,end)));
        xlim([0,0.2]); xlabel('mean speed (um/s)');

        subplot(n_conditions,n_subplot_columns,base+7); hold on;
        histogram(concat_total_distance, 25, 'BinLimits', [0,15], 'Normalization','probability', 'FaceColor', assemble_colors(min(i_condition,end)));
        xlim([0,15]); xlabel('total distance (um)');
    end

    out_png8 = fullfile('graphs', experiment_date, 'assemble_data.png');
    try
        formatfig(f8, 'presentation', out_png8);
    catch
        saveas(f8, out_png8);
    end

    %% Fractional / boxplot summary across conditions
    % Prepare padded matrices for boxplots: fractional long-distance and a few
    % per-file summary metrics (mean maximum displacement, mean total distance, mean speed)
    max_length = max(cellfun(@length, assemble_track_fractional_long_distance));
    
    % If max_length = 1, boxplots wont be meaningful, so we skip
    if max_length > 1
        fractional_distance_matrix = nan(max_length, n_conditions);
        max_displacement_matrix = nan(max_length, n_conditions);
        mean_total_distance_matrix = nan(max_length, n_conditions);
        std_total_distance_matrix = nan(max_length, n_conditions);
        mean_speed_matrix = nan(max_length, n_conditions);

        for i_condition = 1:n_conditions
            % fractional distances is already per-file scalar values
            fractional_vals = assemble_track_fractional_long_distance{i_condition};
            fractional_distance_matrix(1:length(fractional_vals), i_condition) = fractional_vals(:);

            % now build per-file summary stats
            for i_file = 1:length(assemble_track_total_distance{i_condition})
                max_displacement_matrix(i_file,i_condition) = mean(assemble_track_maximum_distance{i_condition}{i_file});
                mean_total_distance_matrix(i_file,i_condition) = mean(assemble_track_total_distance{i_condition}{i_file});
                std_total_distance_matrix(i_file,i_condition) = std(assemble_track_total_distance{i_condition}{i_file});
                mean_speed_matrix(i_file,i_condition) = mean(assemble_track_mean_speed{i_condition}{i_file});
            end
        end

        f9 = figure('Name','Assemble fractional distance','NumberTitle','off'); clf; hold on;
        f9.Position(3:4) = [1200,400];

        subplot(1,5,1);
        boxplot(fractional_distance_matrix, 'Labels', experiment_conditions);
        ylim([0, .3]); ylabel({'fraction of tracks with','maximum track displacement > 1 um'});

        subplot(1,5,2);
        boxplot(max_displacement_matrix, 'Labels', experiment_conditions);
        ylim([0, 1]); ylabel('maximum track displacement (um)');

        subplot(1,5,3);
        boxplot(mean_total_distance_matrix, 'Labels', experiment_conditions);
        ylim([0, 5]); ylabel('mean total track distance (um)');

        subplot(1,5,4);
        boxplot(std_total_distance_matrix, 'Labels', experiment_conditions);
        ylim([0, 5]); ylabel('std total track distance (um)');

        subplot(1,5,5);
        boxplot(mean_speed_matrix, 'Labels', experiment_conditions);
        ylim([0, 0.2]); ylabel('mean speed (um/s)');

        out_png9 = fullfile('graphs', experiment_date, 'assemble_fractional_distance.png');
        try
            formatfig(f9, 'presentation', out_png9);
        catch
            saveas(f9, out_png9);
        end

    else
        warning('Not enough files per condition to create boxplots (max_length = %d). Skipping boxplot summary.', max_length);
    end

    %% Save assembled data to MAT file for downstream analysis
    save(fullfile('graphs', experiment_date, 'assemble_track_data.mat'), ...
        'assemble_track_full', ...
        'assemble_track_lengths', ...
        'assemble_track_velocities', ...
        'assemble_track_mean_speed', ...
        'assemble_track_displacements', ...
        'assemble_track_total_distance', ...
        'assemble_track_maximum_distance', ...
        'assemble_track_mean_squared_displacement', ...
        'assemble_track_fractional_long_distance', ...
        'assemble_track_summary_counts');
end

%% Local helper functions
function ensureDir(folder)
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end