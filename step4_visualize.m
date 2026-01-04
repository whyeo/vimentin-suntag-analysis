% step4_visualize.m
% -------------------------------------------------------------------------
% Title:    step4_visualize
% Purpose:  Create summary visualizations from assembled track data saved by
%           previous steps (assemble_track_data.mat).
%
% Usage:    Edit 'experiment_date' and 'experiment_conditions' below and run
%           this script from the project root.
%
% Inputs:   graphs/<experiment_date>/assemble_track_data.mat
% Outputs:  graphs/<experiment_date>/*.png
%
% Notes:    This script assumes that the data has already been processed by
%           previous steps and saved in the specified location.
%           If there are additional data metrics to visualize, extend the code
%           by uncommenting and modifying the relevant sections below.
%
% Author:   Wei Hong Yeo
% Date:     2026-01-04
% -------------------------------------------------------------------------

% load the saved data file
experiment_date = '2023-12-17';
% list of experimental condition folders under examples/<date>/
experiment_conditions = {'Condition1', 'Condition2'};

load(fullfile('graphs', experiment_date, 'assemble_track_data.mat'),...
    'assemble_track_full',...
    'assemble_track_lengths',...
    'assemble_track_velocities',...
    'assemble_track_mean_speed',...
    'assemble_track_displacements',...
    'assemble_track_total_distance',...
    'assemble_track_maximum_distance',...
    'assemble_track_mean_squared_displacement',...
    'assemble_track_fractional_long_distance',...
    'assemble_track_total_number');

%% for each condition, save the raw assemble data into an excel sheet with
% each sheet corresponding to a file

% create a new assemble data
track_fractional_long_distance = nan(length(experiment_conditions),50);
track_fractional_short_distance = nan(length(experiment_conditions),50);

for i_condition = 1:length(experiment_conditions)
    % get a list of the files within the directory of each condition
    filelist = struct2table(dir(fullfile('examples', experiment_date, experiment_conditions{i_condition}, 'tsoutput', '*.csv')));

    for i_file = 1:length(assemble_track_displacements{i_condition})
        % save the raw data into an excel sheet
        ensureDir(fullfile('graphs', experiment_date));

        track_fractional_short_distance(i_condition,i_file) = sum((assemble_track_displacements{i_condition}{i_file} <= 1) & assemble_track_displacements{i_condition}{i_file} > 0.3)/assemble_track_total_number{i_condition}(i_file);
        track_fractional_long_distance(i_condition,i_file) = sum(assemble_track_displacements{i_condition}{i_file} > 1)/assemble_track_total_number{i_condition}(i_file);

    %     filename = sprintf('graphs/%s/data_track_total_distance.xlsx',experiment_date);
    %     if exist(filename, 'file') == 2
    %         writematrix(assemble_track_total_distance{i_condition}{i_file},filename,...
    %             'Sheet',filelist.name{i_file}(1:end-8),...
    %             'WriteMode','append');
    %     else
    %         writematrix(assemble_track_total_distance{i_condition}{i_file},filename,...
    %             'Sheet',filelist.name{i_file}(1:end-8));
    %     end

    %     filename = sprintf('graphs/%s/data_track_maximum_distance.xlsx',experiment_date);
    %     if exist(filename, 'file') == 2
    %         writematrix(assemble_track_maximum_distance{i_condition}{i_file},filename,...
    %             'Sheet',filelist.name{i_file}(1:end-8),...
    %             'WriteMode','append');
    %     else
    %         writematrix(assemble_track_maximum_distance{i_condition}{i_file},filename,...
    %             'Sheet',filelist.name{i_file}(1:end-8));
    %     end

    %     filename = sprintf('graphs/%s/data_track_mean_speed.xlsx',experiment_date);
    %     if exist(filename, 'file') == 2
    %         writematrix(assemble_track_mean_speed{i_condition}{i_file},filename,...
    %             'Sheet',filelist.name{i_file}(1:end-8),...
    %             'WriteMode','append');
    %     else
    %         writematrix(assemble_track_mean_speed{i_condition}{i_file},filename,...
    %             'Sheet',filelist.name{i_file}(1:end-8));
    %     end

    %     filename = sprintf('graphs/%s/data_track_velocities.xlsx',experiment_date);
    %     if exist(filename, 'file') == 2
    %         writematrix(assemble_track_velocities{i_condition}{i_file},filename,...
    %             'Sheet',filelist.name{i_file}(1:end-8),...
    %             'WriteMode','append');
    %     else
    %         writematrix(assemble_track_velocities{i_condition}{i_file},filename,...
    %             'Sheet',filelist.name{i_file}(1:end-8));
    %     end
    end

    % % save the fractional short distance data
    % filename = sprintf('graphs/%s/data_track_fractional_short_distance.xlsx',experiment_date);
    % if exist(filename, 'file') == 2
    %     writematrix(assemble_track_fractional_short_distance{i_condition},filename,...
    %         'Sheet',experiment_condition(i_condition),...
    %         'WriteMode','append');
    % else
    %     writematrix(assemble_track_fractional_short_distance{i_condition},filename,...
    %         'Sheet',experiment_condition(i_condition));
    % end

    % % save the fractional long distance data
    % filename = sprintf('graphs/%s/data_track_fractional_long_distance.xlsx',experiment_date);
    % if exist(filename, 'file') == 2
    %     writematrix(assemble_track_fractional_long_distance{i_condition},filename,...
    %         'Sheet',experiment_condition(i_condition),...
    %         'WriteMode','append');
    % else
    %     writematrix(assemble_track_fractional_long_distance{i_condition},filename,...
    %         'Sheet',experiment_condition(i_condition));
    % end

    % track_lengths = zeros(sum(cellfun(@length,assemble_track_lengths{i_condition})),1);
    % track_velocities = zeros(sum(cellfun(@length,assemble_track_velocities{i_condition})),1);
    % track_mean_speed = zeros(sum(cellfun(@length,assemble_track_mean_speed{i_condition})),1);
    % track_displacements = zeros(sum(cellfun(@length,assemble_track_displacements{i_condition})),1);
    % track_total_distance = zeros(sum(cellfun(@length,assemble_track_total_distance{i_condition})),1);
    % track_maximum_distance = zeros(sum(cellfun(@length,assemble_track_maximum_distance{i_condition})),1);
    % track_mean_squared_displacement = zeros(sum(cellfun(@length,assemble_track_mean_squared_displacement{i_condition})),1);
    % i_start = 1;

    % for i_file = 1:length(assemble_track_lengths{i_condition})
    %     track_lengths(i_start:i_start+length(assemble_track_lengths{i_condition}{i_file})-1) = assemble_track_lengths{i_condition}{i_file};
    %     track_velocities(i_start:i_start+length(assemble_track_velocities{i_condition}{i_file})-1) = assemble_track_velocities{i_condition}{i_file};
    %     track_mean_speed(i_start:i_start+length(assemble_track_mean_speed{i_condition}{i_file})-1) = assemble_track_mean_speed{i_condition}{i_file};
    %     track_displacements(i_start:i_start+length(assemble_track_displacements{i_condition}{i_file})-1) = assemble_track_displacements{i_condition}{i_file};
    %     track_total_distance(i_start:i_start+length(assemble_track_total_distance{i_condition}{i_file})-1) = assemble_track_total_distance{i_condition}{i_file};
    %     track_maximum_distance(i_start:i_start+length(assemble_track_maximum_distance{i_condition}{i_file})-1) = assemble_track_maximum_distance{i_condition}{i_file};
    %     track_mean_squared_displacement(i_start:i_start+length(assemble_track_mean_squared_displacement{i_condition}{i_file})-1) = assemble_track_mean_squared_displacement{i_condition}{i_file};
    %     i_start = i_start + length(assemble_track_lengths{i_condition}{i_file});
    % end

    % filename = sprintf('graphs/%s/data_track_fractional_long_distance.xlsx',experiment_date);
    % if exist(filename, 'file') == 2
    %     writematrix(assemble_track_fractional_long_distance{i_condition},filename,...
    %         'Sheet',experiment_condition(i_condition),...
    %         'WriteMode','append');
    % else
    %     writematrix(assemble_track_fractional_long_distance{i_condition},filename,...
    %         'Sheet',experiment_condition(i_condition));
    % end

    % filename = sprintf('graphs/%s/data_track_mean_mean_speed.xlsx',experiment_date);
    % if exist(filename, 'file') == 2
    %     writematrix(assemble_track_mean_mean_speed{i_condition},filename,...
    %         'Sheet',experiment_condition(i_condition),...
    %         'WriteMode','append');
    % else
    %     writematrix(assemble_track_mean_mean_speed{i_condition},filename,...
    %         'Sheet',experiment_condition(i_condition));
    % end

    % filename = sprintf('graphs/%s/data_track_mean_total_distance.xlsx',experiment_date);
    % if exist(filename, 'file') == 2
    %     writematrix(assemble_track_mean_total_distance{i_condition},filename,...
    %         'Sheet',experiment_condition(i_condition),...
    %         'WriteMode','append');
    % else
    %     writematrix(assemble_track_mean_total_distance{i_condition},filename,...
    %         'Sheet',experiment_condition(i_condition));
    % end

    % filename = sprintf('graphs/%s/data_track_mean_maximum_distance.xlsx',experiment_date);
    % if exist(filename, 'file') == 2
    %     writematrix(assemble_track_mean_maximum_distance{i_condition},filename,...
    %         'Sheet',experiment_condition(i_condition),...
    %         'WriteMode','append');
    % else
    %     writematrix(assemble_track_mean_maximum_distance{i_condition},filename,...
    %         'Sheet',experiment_condition(i_condition));
    % end
end

% plot a stacked bar graph of the fractional long and short distances

colors_accent = colors('accent');
col_red = colors_accent(1,:);
col_blue = colors_accent(3,:);
col_gray = [0.5,0.5,0.5];

figure(1);
tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');

nexttile;
boxplot(track_fractional_long_distance', 'Labels', experiment_conditions);
title('Fractional Long Distance');
xlabel('Condition');
ylabel('Fractional Long Distance');

nexttile;
boxplot(track_fractional_short_distance', 'Labels', experiment_conditions);
title('Fractional Short Distance');
xlabel('Condition');
ylabel('Fractional Short Distance');

nexttile;
stacked_bar = [1-mean(track_fractional_long_distance,2,'omitnan')-mean(track_fractional_short_distance,2,'omitnan'), mean(track_fractional_short_distance,2,'omitnan'), mean(track_fractional_long_distance,2,'omitnan')];

b = bar(stacked_bar, 'stacked');
b(1).FaceColor = col_gray;
b(2).FaceColor = col_red;
b(3).FaceColor = col_blue;
xlabel('Condition');
xticklabels(experiment_conditions);
ylabel('Fractional Distance');

%% Local helper
function ensureDir(folder)
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end