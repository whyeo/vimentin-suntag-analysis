% step2_analyze_per_file.m
% -------------------------------------------------------------------------
% Title:    step2_analyze_per_file
% Purpose:  Load per-file track caches produced by step1, compute per-track
%           metrics (length, displacement, total distance, mean speed, MSD,
%           directional persistence, alpha exponent) and produce per-file
%           diagnostic figures. Save a per-file summary .mat for aggregation.
%
% Usage:    Edit the 'experiment_date' and 'experiment_conditions' variables
%           below and run this script from the project root.
%
% Inputs:   examples/<experiment_date>/<condition>/tsoutput/*_track.mat (caches)
% Outputs:  graphs/<experiment_date>/<condition>/*_MSD_DP.mat, PNG figures
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

% color helper (optional)
try
    colors_accent = colors('accent');
catch
    colors_accent = [0.2,0.6,0.9; 0.9,0.35,0.2; 0.6,0.3,0.9]; % fallback colors
end

% track type labels
STATIONARY_TRACK = 0;
SHORT_TRACK = 1;
LONG_TRACK = 2;

% initialize containers (one cell per condition)
n_conditions = numel(experiment_conditions);
assemble_track_full = cell(1, n_conditions);
assemble_track_lengths = cell(1, n_conditions);
assemble_track_velocities = cell(1, n_conditions);
assemble_track_displacements = cell(1, n_conditions);
assemble_track_mean_speed = cell(1, n_conditions);
assemble_track_total_distance = cell(1, n_conditions);
assemble_track_maximum_distance = cell(1, n_conditions);
assemble_track_mean_squared_displacement = cell(1, n_conditions);
assemble_track_directional_ratio = cell(1, n_conditions);
assemble_track_directional_autocorrelation = cell(1, n_conditions);

assemble_track_fractional_long_distance = cell(1, n_conditions);
assemble_track_fractional_short_distance = cell(1, n_conditions);
assemble_track_fractional_stationary = cell(1, n_conditions);
assemble_track_mean_maximum_distance = cell(1, n_conditions);
assemble_track_mean_total_distance = cell(1, n_conditions);
assemble_track_mean_mean_speed = cell(1, n_conditions);
assemble_track_total_number = cell(1, n_conditions);

% ensure top-level output folder exists
ensureDir(fullfile('graphs', experiment_date));

%% Process each condition and each file
for i_condition = 1:n_conditions
    cond_name = experiment_conditions{i_condition};
    ts_folder = fullfile('examples', experiment_date, cond_name, 'tsoutput');

    % find ThunderSTORM CSVs / track caches for this condition
    filelist = dir(fullfile(ts_folder, '*-sub.csv'));
    n_files = numel(filelist);

    % preallocate per-condition storage
    assemble_track_full{i_condition} = cell(1, n_files);
    assemble_track_lengths{i_condition} = cell(1, n_files);
    assemble_track_velocities{i_condition} = cell(1, n_files);
    assemble_track_mean_speed{i_condition} = cell(1, n_files);
    assemble_track_displacements{i_condition} = cell(1, n_files);
    assemble_track_total_distance{i_condition} = cell(1, n_files);
    assemble_track_maximum_distance{i_condition} = cell(1, n_files);
    assemble_track_mean_squared_displacement{i_condition} = cell(1, n_files);
    assemble_track_directional_ratio{i_condition} = cell(1, n_files);
    assemble_track_directional_autocorrelation{i_condition} = cell(1, n_files);

    assemble_track_fractional_long_distance{i_condition} = zeros(1, n_files);
    assemble_track_fractional_short_distance{i_condition} = zeros(1, n_files);
    assemble_track_fractional_stationary{i_condition} = zeros(1, n_files);
    assemble_track_mean_maximum_distance{i_condition} = zeros(1, n_files);
    assemble_track_mean_total_distance{i_condition} = zeros(1, n_files);
    assemble_track_mean_mean_speed{i_condition} = zeros(1, n_files);
    assemble_track_total_number{i_condition} = zeros(1, n_files);

    % ensure per-condition output folder
    cond_out = fullfile('graphs', experiment_date, cond_name);
    ensureDir(cond_out);

    for i_file = 1:n_files
        csv_name = filelist(i_file).name;
        csv_path = fullfile(ts_folder, csv_name);
        fprintf('Processing %s / %s\n', cond_name, csv_name);

        % corresponding track cache (produced by step1_tracking)
        track_cache = fullfile(ts_folder, sprintf('%s_track.mat', csv_name(1:end-8)));
        if ~exist(track_cache, 'file')
            warning('  Track cache not found: %s — skipping file', track_cache);
            continue;
        end

        % load tracks variable
        S = load(track_cache, 'tracks');
        tracks = S.tracks; % columns expected: x, y, intensity, frame, ..., track_id

        % basic bookkeeping
        track_max_length = max(tracks(:,end-1));
        n_tracks = max(tracks(:,end));

        % per-file containers
        track_full = cell(n_tracks,1);
        track_length = zeros(n_tracks,1);
        track_velocities = nan(n_tracks,1);
        track_mean_speed = nan(n_tracks,1);
        track_displacement = zeros(n_tracks,1);
        track_total_distance = zeros(n_tracks,1);
        track_maximum_distance = zeros(n_tracks,1);
        track_maximum_velocity = zeros(n_tracks,1);
        track_type = zeros(n_tracks,1);

        track_mean_squared_displacement = nan(n_tracks, track_max_length);
        track_alpha = nan(n_tracks,1);
        track_directional_ratios = nan(n_tracks, track_max_length);
        track_directional_autocorrelation = nan(n_tracks, track_max_length);

        % iterate over tracks
        for i_track = 1:n_tracks
            curr = tracks(tracks(:,end) == i_track, :);
            track_full{i_track} = curr;

            % frame-based length (number of frames spanned)
            frames = curr(:, end-1);
            length_frames = max(frames) - min(frames);
            track_length(i_track) = length_frames * frame_interval_s;

            % displacement and distances in pixels
            displacement_px = norm(curr(end,1:2) - curr(1,1:2));
            total_dist_px = sum(sqrt(sum(diff(curr(:,1:2)).^2,2)));
            max_step_px = max(sqrt(sum(diff(curr(:,1:2)).^2,2)));

            % convert to physical units (um)
            displacement_um = displacement_px * (pixel_size_nm/1000);
            total_dist_um = total_dist_px * (pixel_size_nm/1000);
            max_step_um = max_step_px * (pixel_size_nm/1000);

            % velocities (guard against zero-length tracks)
            if length_frames > 0
                velocity_um_s = (displacement_um / (length_frames * frame_interval_s));
                mean_speed_um_s = (total_dist_um / (length_frames * frame_interval_s));
            else
                velocity_um_s = NaN;
                mean_speed_um_s = NaN;
            end

            % store
            track_velocities(i_track) = velocity_um_s;
            track_mean_speed(i_track) = mean_speed_um_s;
            track_displacement(i_track) = displacement_um;
            track_total_distance(i_track) = total_dist_um;
            track_maximum_distance(i_track) = max_step_um;
            track_maximum_velocity(i_track) = max_step_um / frame_interval_s;

            % classify
            if track_maximum_distance(i_track) > 1
                track_type(i_track) = LONG_TRACK;
            elseif track_maximum_distance(i_track) > 0.3
                track_type(i_track) = SHORT_TRACK;
            else
                track_type(i_track) = STATIONARY_TRACK;
            end

            % prepare pairwise distances and dot-products for MSD and directionality
            n_points = length_frames + 1; % inclusive
            % build an xy matrix aligned to contiguous frames (missing frames padded with NaN)
            xy = nan(n_points, 2);
            start_frame = min(frames);
            for p = 1:n_points
                fr = start_frame + p - 1;
                idx = frames == fr;
                if any(idx)
                    xy(p,:) = curr(idx, 1:2);
                end
            end

            % distances matrix (upper triangle)
            dx = triu(xy(:,1) - xy(:,1)');
            dy = triu(xy(:,2) - xy(:,2)');
            distances = sqrt(dx.^2 + dy.^2);

            % directional dot-products for autocorrelation (sparse estimate)
            v = [diff(xy(:,1)), diff(xy(:,2))];
            % compute dot-products between successive steps where available
            if size(v,1) >= 1
                step_dots = NaN(size(v,1), size(v,1));
                for a = 1:size(v,1)
                    for b = 1:size(v,1)
                        sd = dot(v(a,:), v(b,:));
                        denom = norm(v(a,:)) * norm(v(b,:));
                        if denom ~= 0
                            step_dots(a,b) = sd / denom;
                        else
                            step_dots(a,b) = NaN;
                        end
                    end
                end
            else
                step_dots = NaN(1,1);
            end

            % compute MSD and directional metrics across lags
            for lag = 1:max(1, n_points-1)
                % mean squared displacement for this lag (ignore NaNs)
                diag_vals = diag(distances, lag);
                if ~isempty(diag_vals)
                    track_mean_squared_displacement(i_track, lag) = mean(diag_vals, 'omitnan');
                else
                    track_mean_squared_displacement(i_track, lag) = NaN;
                end

                % directional ratio d/D (simple definition)
                if lag == 1
                    track_directional_ratios(i_track, lag) = 1;
                else
                    d = distances(1, lag+1);
                    D = sum(diag(distances(1:lag+1, 1:lag+1), 1), 'omitnan');
                    if D ~= 0
                        track_directional_ratios(i_track, lag) = d / D;
                    else
                        track_directional_ratios(i_track, lag) = NaN;
                    end
                end

                % directional autocorrelation (mean of step dot-products at given lag)
                if lag <= size(step_dots,1)
                    % use mean of off-diagonal at distance 'lag'
                    track_directional_autocorrelation(i_track, lag) = mean(diag(step_dots, lag), 'omitnan');
                else
                    track_directional_autocorrelation(i_track, lag) = NaN;
                end
            end

            % compute alpha from MSD if enough points
            msd_vals = track_mean_squared_displacement(i_track, :);
            valid = ~isnan(msd_vals);
            if sum(valid) > 3
                lags = (1:sum(valid)) * frame_interval_s;
                try
                    p = polyfit(log(lags), log(msd_vals(valid)), 1);
                    track_alpha(i_track) = p(1);
                catch
                    track_alpha(i_track) = NaN;
                end
            end
        end % per-track

        % per-file summary fractions
        track_fractional_long_distance = sum(track_type == LONG_TRACK) / max(1, n_tracks);
        track_fractional_short_distance = sum(track_type == SHORT_TRACK) / max(1, n_tracks);
        track_fractional_stationary = sum(track_type == STATIONARY_TRACK) / max(1, n_tracks);

        % convert MSD to um^2 (already in pixel units)
        track_mean_squared_displacement = track_mean_squared_displacement * (pixel_size_nm/1000)^2;
        timelags = (0:track_max_length-1) * frame_interval_s;

        % optional plotting (MSD, directional persistence, distributions)
        try
            f1 = figure('Name', sprintf('%s - %s', cond_name, csv_name(1:end-8)), 'NumberTitle','off');
            f1.Position(3:4) = [1600, 900];
            tiledlayout(5,3,'TileSpacing','Compact','Padding','Compact');

            % panel: tracks
            nexttile(1,[5,1]); hold on; box on; axis equal;
            for it = 1:n_tracks
                col = [0.5,0.5,0.5];
                if track_type(it) == LONG_TRACK, col = colors_accent(3,:); end
                if track_type(it) == SHORT_TRACK, col = colors_accent(1,:); end
                plot(track_full{it}(:,1), track_full{it}(:,2), 'Color', [col, 0.5]);
            end

            % panel: raw MSD traces (log-log)
            nexttile(2); hold on; box on; set(gca,'XScale','log','YScale','log'); xlabel('Time lag (s)'); ylabel('MSD (\mum^2)');
            for it = 1:n_tracks
                plot(timelags(1:size(track_mean_squared_displacement,2)), track_mean_squared_displacement(it, :), 'Color', [0.3,0.3,0.3,0.3]);
            end

            % panel: mean MSD +/- std per class
            nexttile(3); hold on; box on; xlabel('Time lag (s)'); ylabel('MSD (\mum^2)');
            errorbar(timelags, mean(track_mean_squared_displacement(track_type == LONG_TRACK, :), 1, 'omitnan'), std(track_mean_squared_displacement(track_type == LONG_TRACK, :), 0, 1, 'omitnan'), 'Color', colors_accent(3,:));
            errorbar(timelags, mean(track_mean_squared_displacement(track_type == SHORT_TRACK, :), 1, 'omitnan'), std(track_mean_squared_displacement(track_type == SHORT_TRACK, :), 0, 1, 'omitnan'), 'Color', colors_accent(1,:));
            errorbar(timelags, mean(track_mean_squared_displacement(track_type == STATIONARY_TRACK, :), 1, 'omitnan'), std(track_mean_squared_displacement(track_type == STATIONARY_TRACK, :), 0, 1, 'omitnan'), 'Color', [0.5,0.5,0.5]);

            % panel: directional persistence
            nexttile(5); hold on; box on; xlabel('Time lag (s)'); ylabel('Directional Persistence (d/D)');
            for it = 1:n_tracks
                plot(timelags(1:size(track_directional_ratios,2)), track_directional_ratios(it, :), 'Color', [0.6,0.6,0.6,0.3]);
            end

            % panel: mean directional persistence per class
            nexttile(6); hold on; box on; xlabel('Time lag (s)'); ylabel('Directional Persistence (d/D)');
            errorbar(timelags, mean(track_directional_ratios(track_type == LONG_TRACK, :), 1, 'omitnan'), std(track_directional_ratios(track_type == LONG_TRACK, :), 0, 1, 'omitnan'), 'Color', colors_accent(3,:));
            errorbar(timelags, mean(track_directional_ratios(track_type == SHORT_TRACK, :), 1, 'omitnan'), std(track_directional_ratios(track_type == SHORT_TRACK, :), 0, 1, 'omitnan'), 'Color', colors_accent(1,:));
            errorbar(timelags, mean(track_directional_ratios(track_type == STATIONARY_TRACK, :), 1, 'omitnan'), std(track_directional_ratios(track_type == STATIONARY_TRACK, :), 0, 1, 'omitnan'), 'Color', [0.5,0.5,0.5]);

            % panel: alpha histogram
            nexttile(8); hold on; box on; xlabel('Alpha');
            histogram(track_alpha(track_type == LONG_TRACK), -1:0.05:1.5, 'Normalization', 'probability', 'FaceColor', colors_accent(3,:));
            histogram(track_alpha(track_type == SHORT_TRACK), -1:0.05:1.5, 'Normalization', 'probability', 'FaceColor', colors_accent(1,:));
            histogram(track_alpha(track_type == STATIONARY_TRACK), -1:0.05:1.5, 'Normalization', 'probability', 'FaceColor', [0.5,0.5,0.5]);

            % panel: directional autocorr at a short lag (example)
            nexttile(9); hold on; box on; xlabel('Directional Autocorrelation (cos(\theta)) at lag 3');
            sel_lag = min(3, size(track_directional_autocorrelation,2));
            histogram(track_directional_autocorrelation(track_type == LONG_TRACK, sel_lag), -1:0.05:1, 'Normalization', 'probability', 'FaceColor', colors_accent(3,:));
            histogram(track_directional_autocorrelation(track_type == SHORT_TRACK, sel_lag), -1:0.05:1, 'Normalization', 'probability', 'FaceColor', colors_accent(1,:));
            histogram(track_directional_autocorrelation(track_type == STATIONARY_TRACK, sel_lag), -1:0.05:1, 'Normalization', 'probability', 'FaceColor', [0.5,0.5,0.5]);

            % panel: velocity distributions
            nexttile(14); hold on; box on; xlabel('Mean speed (\mum/s)');
            histogram(track_mean_speed(track_type == LONG_TRACK), 0:0.01:0.6, 'Normalization', 'probability', 'FaceColor', colors_accent(3,:));
            histogram(track_mean_speed(track_type == SHORT_TRACK), 0:0.01:0.6, 'Normalization', 'probability', 'FaceColor', colors_accent(1,:));
            histogram(track_mean_speed(track_type == STATIONARY_TRACK), 0:0.01:0.6, 'Normalization', 'probability', 'FaceColor', [0.5,0.5,0.5]);

            % panel: max speed distributions
            nexttile(15); hold on; box on; xlabel('Maximum speed (\mum/s)');
            histogram(track_maximum_velocity(track_type == LONG_TRACK), 0:0.02:1.2, 'Normalization', 'probability', 'FaceColor', colors_accent(3,:));
            histogram(track_maximum_velocity(track_type == SHORT_TRACK), 0:0.02:1.2, 'Normalization', 'probability', 'FaceColor', colors_accent(1,:));
            histogram(track_maximum_velocity(track_type == STATIONARY_TRACK), 0:0.02:1.2, 'Normalization', 'probability', 'FaceColor', [0.5,0.5,0.5]);

            % save figure (try formatfig, else fallback)
            out_png = fullfile(cond_out, sprintf('%s_MSD_DP.png', csv_name(1:end-8)));
            try
                formatfig(f1, 'publication', out_png);
            catch
                saveas(f1, out_png);
            end
            close(f1);
        catch ME
            warning('  Could not create figure for %s: %s', csv_name, ME.message);
        end

        % Save per-file summary .mat
        save_data_file = fullfile(cond_out, sprintf('%s_MSD_DP.mat', csv_name(1:end-8)));
        try
            save(save_data_file, 'tracks', 'track_mean_squared_displacement', 'track_alpha', 'track_directional_ratios', 'track_directional_autocorrelation', 'timelags', 'track_type', 'track_fractional_long_distance', 'track_fractional_short_distance', 'track_fractional_stationary','track_mean_speed', 'track_total_distance', 'track_maximum_distance', 'track_maximum_velocity');
        catch ME
            warning('  Could not save MAT for %s: %s', csv_name, ME.message);
        end

        % Store assembled per-file results for later aggregation
        assemble_track_full{i_condition}{i_file} = track_full;
        assemble_track_lengths{i_condition}{i_file} = track_length;
        assemble_track_velocities{i_condition}{i_file} = track_velocities;
        assemble_track_mean_speed{i_condition}{i_file} = track_mean_speed;
        assemble_track_displacements{i_condition}{i_file} = track_displacement;
        assemble_track_total_distance{i_condition}{i_file} = track_total_distance;
        assemble_track_maximum_distance{i_condition}{i_file} = track_maximum_distance;
        assemble_track_mean_squared_displacement{i_condition}{i_file} = track_mean_squared_displacement;

        assemble_track_fractional_long_distance{i_condition}(i_file) = track_fractional_long_distance;
        assemble_track_fractional_short_distance{i_condition}(i_file) = track_fractional_short_distance;
        assemble_track_fractional_stationary{i_condition}(i_file) = track_fractional_stationary;
        assemble_track_mean_maximum_distance{i_condition}(i_file) = mean(track_maximum_distance, 'omitnan');
        assemble_track_mean_total_distance{i_condition}(i_file) = mean(track_total_distance, 'omitnan');
        assemble_track_mean_mean_speed{i_condition}(i_file) = mean(track_mean_speed, 'omitnan');
        assemble_track_total_number{i_condition}(i_file) = n_tracks;

    end % per-file
end % per-condition

%% Local helper
function ensureDir(folder)
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end