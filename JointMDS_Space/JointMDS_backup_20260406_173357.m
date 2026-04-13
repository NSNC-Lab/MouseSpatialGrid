clearvars
close all

% Load in the data and simulation.
analysis_dir = 'C:\Users\ipboy\Desktop\SingleChannelFigures\JointMDS_Space';
cd(analysis_dir)

load('50epoch_100ms_and_300ms_loss_All_cells.mat')
load('all_units_info_with_polished_criteria_modified_perf.mat', 'all_data');

% Analysis settings.
bin_width = 200;
sample_rate_hz = 10000;
time_len = 29801;
n_cells = size(output, 1);
n_trials = size(output, 3);
n_bins = floor(time_len / bin_width);
zoom_xlim = [-35, -15];
zoom_ylim = [-10, 10];
stimulus_path = 'C:\Users\ipboy\Documents\GitHub\MouseSpatialGrid\LearningModels\resampled-stimuli\target\200k_target1.wav';

feature_display_names = { ...
    'Mean Rate (Hz)'; ...
    'Peak Rate (1st on/off, Hz)'; ...
    'Peak Latency (1st on/off, ms)'; ...
    'Onset Latency (1st on/off, 20% peak, ms)'; ...
    'Half-Max Width (1st on/off, ms)'; ...
    'Temporal Center (ms)'; ...
    'log2(Early/Late)'; ...
    'Temporal Sparsity'; ...
    'Trial Reliability'; ...
    'Peak-Bin Fano (1st on/off)' ...
    };

[target1, stimulus_fs] = audioread(stimulus_path);
stimulus_event_info = get_stimulus_event_info(target1, stimulus_fs, sample_rate_hz, bin_width, n_bins);
plot_stimulus_event_windows(stimulus_event_info, 'Target1_first_onset_offset_windows.pdf')

% Extract out data and simulation rasters/PSTHs.
Rasters_sim = zeros(n_cells, n_trials, time_len);
Rasters_data = zeros(n_cells, n_trials, time_len);
PSTHs_sim = zeros(n_cells, n_bins);
PSTHs_data = zeros(n_cells, n_bins);
features_sim = nan(n_cells, numel(feature_display_names));
features_data = nan(n_cells, numel(feature_display_names));

final_binned_sse = squeeze(losses(end, 2, :, :)).';
[~, idx] = min(final_binned_sse, [], 1);

for k = 1:n_cells
    sim_raster = squeeze(output(k, idx(k), :, :, :));
    sim_raster = squeeze(sim_raster);
    if isvector(sim_raster)
        sim_raster = reshape(sim_raster, n_trials, []);
    elseif size(sim_raster, 1) ~= n_trials && size(sim_raster, 2) == n_trials
        sim_raster = sim_raster.';
    elseif size(sim_raster, 1) ~= n_trials
        sim_raster = reshape(sim_raster, n_trials, []);
    end
    sim_raster = double(sim_raster);
    Rasters_sim(k, :, :) = sim_raster;

    sim_trial_bins = raster_to_trial_bins(sim_raster, bin_width);
    PSTHs_sim(k, :) = sum(sim_trial_bins, 1);
    features_sim(k, :) = compute_response_features(sim_trial_bins, sample_rate_hz, bin_width, stimulus_event_info);

    focus = get_focus_index(all_data(k).tuning_type);
    spike_times = all_data(k).ctrl_tar1_timestamps(:, focus);
    data_raster = zeros(n_trials, time_len);
    for m = 1:n_trials
        stim_mask = (spike_times{m} > 0) & (spike_times{m} < 2.9801);
        trial_indices = round(spike_times{m}(stim_mask) * sample_rate_hz);
        trial_indices = trial_indices(trial_indices >= 0 & trial_indices < time_len);
        data_raster(m, trial_indices + 1) = 1;
    end
    Rasters_data(k, :, :) = data_raster;

    data_trial_bins = raster_to_trial_bins(data_raster, bin_width);
    PSTHs_data(k, :) = sum(data_trial_bins, 1);
    features_data(k, :) = compute_response_features(data_trial_bins, sample_rate_hz, bin_width, stimulus_event_info);
end

% Joint MDS analysis.
PSTHs_all = [PSTHs_data; PSTHs_sim];
comp_matrix = zeros(2 * n_cells, 2 * n_cells);
for m = 1:(2 * n_cells)
    diffs = PSTHs_all - PSTHs_all(m, :);
    comp_matrix(m, :) = sqrt(sum(diffs .^ 2, 2)).';
end

cells = [1:n_cells, 1:n_cells];
[Y, stress] = mdscale(comp_matrix, 2, 'criterion', 'metricsstress');
disp("MDS stress: " + string(stress))

% MDS Plot (1) - Full plot, uncolored.
f = figure('Position', [50, 50, 1000, 1000]);
scatter(Y(:, 1), Y(:, 2), 36, 'bo', 'filled'); hold on
labels = string(cells);
text(Y(:, 1), Y(:, 2), labels, 'FontSize', 10)
title('Full plot - Uncolored')
xlabel('MDS-Axis 1 (Firing Rate)')
ylabel('MDS-Axis 2 (Not Fully known)')

% MDS Plot (2) - Zoomed.
f = figure('Position', [50, 50, 1000, 1000]);
scatter(Y(:, 1), Y(:, 2), 36, 'bo', 'filled'); hold on
text(Y(:, 1), Y(:, 2), labels, 'FontSize', 10)
title('Full plot - Uncolored')
xlabel('MDS-Axis 1 (Firing Rate)')
ylabel('MDS-Axis 2 (Not Fully known)')
ylim([-5, 5])
xlim([-30, 10])

% Color by self-similarity (normalized inverse distance between matching
% data/model points).
distances = zeros(1, n_cells);
for k = 1:n_cells
    dist_vals = 1 ./ comp_matrix(n_cells + k, 1:n_cells);
    dist_vals_norm = dist_vals ./ max(dist_vals);
    distances(k) = dist_vals_norm(k);
end
distances = [distances, distances];

% MDS Plot (3) - Colored by self-similarity.
f = figure('Position', [50, 50, 1000, 1000]);
scatter(Y(:, 1), Y(:, 2), 36, distances, 'filled'); hold on
text(Y(:, 1), Y(:, 2), labels, 'FontSize', 10)
title('Full plot - Colored by distance Matrix')
xlabel('MDS-Axis 1 (Firing Rate)')
ylabel('MDS-Axis 2 (Not Fully known)')
colorbar;

% MDS Plot (4) - Colored with connecting lines.
f = figure('Position', [50, 50, 1000, 1000]);
scatter(Y(:, 1), Y(:, 2), 36, distances, 'filled');
hold on
cmap = parula(256);
line_distances = distances(1:n_cells);
cmin = min(line_distances);
cmax = max(line_distances);
clim([cmin, cmax])

for k = 1:n_cells
    t = (line_distances(k) - cmin) / (cmax - cmin + eps);
    color_idx = max(1, min(size(cmap, 1), round(1 + t * (size(cmap, 1) - 1))));
    this_color = cmap(color_idx, :);

    plot([Y(k, 1), Y(k + n_cells, 1)], ...
         [Y(k, 2), Y(k + n_cells, 2)], ...
         'Color', this_color, 'LineWidth', 1.5);
end

text(Y(:, 1), Y(:, 2), labels, 'FontSize', 10)
title('Full plot - Colored by distance Matrix')
xlabel('MDS-Axis 1 (Firing Rate)')
ylabel('MDS-Axis 2 (Not Fully known)')
colormap(cmap)
colorbar
exportgraphics(f, 'Full_Plot_Colores.pdf', 'ContentType', 'vector')

% MDS Plot (5) - Zoomed, colored with connecting lines.
f = figure('Position', [50, 50, 1000, 1000]);
scatter(Y(:, 1), Y(:, 2), 36, distances, 'filled');
hold on
clim([cmin, cmax])

for k = 1:n_cells
    t = (line_distances(k) - cmin) / (cmax - cmin + eps);
    color_idx = max(1, min(size(cmap, 1), round(1 + t * (size(cmap, 1) - 1))));
    this_color = cmap(color_idx, :);

    plot([Y(k, 1), Y(k + n_cells, 1)], ...
         [Y(k, 2), Y(k + n_cells, 2)], ...
         'Color', this_color, 'LineWidth', 1.5);
end

title('Full plot - Colored by distance Matrix')
xlabel('MDS-Axis 1 (Firing Rate)')
ylabel('MDS-Axis 2 (Not Fully known)')
colormap(cmap)
colorbar
xlim(zoom_xlim)
ylim(zoom_ylim)

mask = Y(:, 1) >= zoom_xlim(1) & Y(:, 1) <= zoom_xlim(2) & ...
       Y(:, 2) >= zoom_ylim(1) & Y(:, 2) <= zoom_ylim(2);
text(Y(mask, 1), Y(mask, 2), labels(mask), 'FontSize', 10, 'Clipping', 'on')
exportgraphics(f, 'Full_Plot_Colored_Zoomed.pdf', 'ContentType', 'vector')

% Feature overlays for explaining the MDS axes.
features_all = [features_data; features_sim];
axis1_corrs_all = nan(numel(feature_display_names), 1);
axis2_corrs_all = nan(numel(feature_display_names), 1);
axis2_corrs_data = nan(numel(feature_display_names), 1);
axis2_corrs_model = nan(numel(feature_display_names), 1);

for f_idx = 1:numel(feature_display_names)
    axis1_corrs_all(f_idx) = safe_corr(Y(:, 1), features_all(:, f_idx));
    axis2_corrs_all(f_idx) = safe_corr(Y(:, 2), features_all(:, f_idx));
    axis2_corrs_data(f_idx) = safe_corr(Y(1:n_cells, 2), features_data(:, f_idx));
    axis2_corrs_model(f_idx) = safe_corr(Y((n_cells + 1):end, 2), features_sim(:, f_idx));
end

feature_summary = table( ...
    string(feature_display_names), ...
    axis1_corrs_all, ...
    axis2_corrs_all, ...
    axis2_corrs_data, ...
    axis2_corrs_model, ...
    'VariableNames', {'Feature', 'Axis1_r_all', 'Axis2_r_all', 'Axis2_r_data', 'Axis2_r_model'});

disp(feature_summary)
writetable(feature_summary, 'JointMDS_feature_axis_correlations.csv')
save('JointMDS_feature_summary.mat', ...
     'Y', 'stress', 'feature_display_names', 'features_data', 'features_sim', ...
     'axis1_corrs_all', 'axis2_corrs_all', 'axis2_corrs_data', 'axis2_corrs_model')

plot_feature_grid(Y, features_all, feature_display_names, axis1_corrs_all, axis2_corrs_all, ...
                  n_cells, [], [], 'JointMDS_feature_grid.pdf')

plot_feature_grid(Y, features_all, feature_display_names, axis1_corrs_all, axis2_corrs_all, ...
                  n_cells, zoom_xlim, zoom_ylim, 'JointMDS_feature_grid_zoomed.pdf')

%% Find verticality of the lines.

verticality = [];
for k = 1:n_cells
    verticality = [verticality, (180 / pi) * atan(abs(Y(k, 2) - Y(k + n_cells, 2)) / abs(Y(k, 1) - Y(k + n_cells, 1)))];
end

figure;
histogram(verticality)
disp('mean angle of lines: ' + string(mean(verticality)))

verticality = [];
boundary = 40;

for k = 1:n_cells
    if hypot(Y(k, 1), Y(k, 2)) > boundary
        if hypot(Y(k + n_cells, 1), Y(k + n_cells, 2)) > boundary
            verticality = [verticality, (180 / pi) * atan(abs(Y(k, 2) - Y(k + n_cells, 2)) / abs(Y(k, 1) - Y(k + n_cells, 1)))];
        end
    end
end

figure;
histogram(verticality)
disp('mean angle of lines: ' + string(mean(verticality)))

function trial_bins = raster_to_trial_bins(raster, bin_width)
    raster = double(raster);
    n_trials_local = size(raster, 1);
    n_bins_local = floor(size(raster, 2) / bin_width);
    trimmed = raster(:, 1:(n_bins_local * bin_width));
    reshaped = reshape(trimmed.', bin_width, n_bins_local, n_trials_local);
    trial_bins = squeeze(sum(reshaped, 1)).';
    if n_trials_local == 1
        trial_bins = reshape(trial_bins, 1, []);
    end
end

function feature_vector = compute_response_features(trial_bins, sample_rate_hz, bin_width, stimulus_event_info)
    bin_sec = bin_width / sample_rate_hz;
    bin_ms = 1000 * bin_sec;
    mean_counts = mean(trial_bins, 1);
    mean_rate = mean_counts / bin_sec;
    n_bins_local = numel(mean_rate);
    bin_centers_ms = ((1:n_bins_local) - 0.5) * bin_ms;

    mean_rate_hz = mean(mean_rate);
    onset_window = clip_indices(stimulus_event_info.onset_bins, n_bins_local);
    offset_window = clip_indices(stimulus_event_info.offset_bins, n_bins_local);
    [peak_rate_hz, peak_idx, chosen_window] = first_event_peak(mean_rate, onset_window, offset_window);

    if peak_rate_hz <= 0 || isempty(chosen_window)
        peak_latency_ms = NaN;
        onset_latency_ms = NaN;
        half_max_width_ms = NaN;
        peak_bin_fano = NaN;
    else
        peak_latency_ms = bin_centers_ms(peak_idx);

        chosen_rate = mean_rate(chosen_window);
        onset_idx_local = find(chosen_rate >= 0.20 * peak_rate_hz, 1, 'first');
        onset_idx = chosen_window(onset_idx_local);
        onset_latency_ms = bin_centers_ms(onset_idx);

        peak_pos_in_window = find(chosen_window == peak_idx, 1, 'first');
        half_max_mask = chosen_rate >= 0.50 * peak_rate_hz;
        left_pos = peak_pos_in_window;
        right_pos = peak_pos_in_window;
        while left_pos > 1 && half_max_mask(left_pos - 1)
            left_pos = left_pos - 1;
        end
        while right_pos < numel(chosen_window) && half_max_mask(right_pos + 1)
            right_pos = right_pos + 1;
        end
        half_max_width_ms = (right_pos - left_pos + 1) * bin_ms;

        peak_bin_counts = trial_bins(:, peak_idx);
        if mean(peak_bin_counts) > 0
            peak_bin_fano = var(peak_bin_counts, 1) / mean(peak_bin_counts);
        else
            peak_bin_fano = NaN;
        end
    end

    total_rate = sum(mean_rate);
    if total_rate > 0
        temporal_center_ms = sum(bin_centers_ms .* mean_rate) / total_rate;
    else
        temporal_center_ms = NaN;
    end

    split_idx = max(1, floor(n_bins_local / 2));
    early_mean = mean(mean_rate(1:split_idx));
    second_half = mean_rate((split_idx + 1):end);
    if isempty(second_half)
        second_half = mean_rate(split_idx);
    end
    late_mean = mean(second_half);
    log2_early_late = log2((early_mean + eps) / (late_mean + eps));

    if all(mean_rate == 0)
        temporal_sparsity = NaN;
    else
        numerator = 1 - ((mean(mean_rate) ^ 2) / (mean(mean_rate .^ 2) + eps));
        denominator = 1 - (1 / n_bins_local);
        temporal_sparsity = numerator / (denominator + eps);
    end

    trial_reliability = mean_pairwise_trial_corr(trial_bins);

    feature_vector = [ ...
        mean_rate_hz, ...
        peak_rate_hz, ...
        peak_latency_ms, ...
        onset_latency_ms, ...
        half_max_width_ms, ...
        temporal_center_ms, ...
        log2_early_late, ...
        temporal_sparsity, ...
        trial_reliability, ...
        peak_bin_fano ...
        ];
end

function stimulus_event_info = get_stimulus_event_info(stimulus_wave, stimulus_fs, response_fs, bin_width, n_bins)
    if size(stimulus_wave, 2) > 1
        stimulus_wave = mean(stimulus_wave, 2);
    end

    stimulus_wave = double(stimulus_wave(:));
    envelope = abs(stimulus_wave);
    smooth_window = max(1, round(0.005 * stimulus_fs));
    envelope = movmean(envelope, smooth_window);

    samples_per_bin = max(1, round(bin_width * stimulus_fs / response_fs));
    needed_samples = n_bins * samples_per_bin;
    if numel(envelope) < needed_samples
        envelope(end + 1:needed_samples) = 0;
    end

    envelope = envelope(1:needed_samples);
    envelope_bins = mean(reshape(envelope, samples_per_bin, n_bins), 1);
    envelope_bins = movmean(envelope_bins, 3);

    if max(envelope_bins) > 0
        envelope_bins = envelope_bins / max(envelope_bins);
    end

    threshold = 0.05;
    start_idx = find(envelope_bins > threshold, 1, 'first');
    if isempty(start_idx)
        start_idx = 1;
    end

    first_peak_idx = find_first_local_max(envelope_bins, start_idx);
    if first_peak_idx <= start_idx
        first_peak_idx = min(n_bins, start_idx + 1);
    end

         = find_first_local_min(envelope_bins, first_peak_idx);
    if first_trough_idx <= first_peak_idx
        first_trough_idx = min(n_bins, first_peak_idx + 1);
    end

    bin_ms = 1000 * (bin_width / response_fs);
    stimulus_event_info = struct();
    stimulus_event_info.envelope_bins = envelope_bins;
    stimulus_event_info.bin_centers_ms = ((1:n_bins) - 0.5) * bin_ms;
    stimulus_event_info.start_idx = start_idx;
    stimulus_event_info.first_peak_idx = first_peak_idx;
    stimulus_event_info.first_trough_idx = first_trough_idx;
    stimulus_event_info.onset_bins = start_idx:first_peak_idx;
    stimulus_event_info.offset_bins = first_peak_idx:first_trough_idx;
end

function peak_idx = find_first_local_max(values, start_idx)
    peak_idx = start_idx;
    for idx = (start_idx + 1):(numel(values) - 1)
        if values(idx) >= values(idx - 1) && values(idx) >= values(idx + 1) && ...
                any(values(start_idx:idx) > values(start_idx))
            peak_idx = idx;
            return
        end
    end

    [~, local_idx] = max(values(start_idx:end));
    peak_idx = start_idx + local_idx - 1;
end

function trough_idx = find_first_local_min(values, start_idx)
    trough_idx = start_idx;
    for idx = (start_idx + 1):(numel(values) - 1)
        if values(idx) <= values(idx - 1) && values(idx) <= values(idx + 1) && ...
                any(values(start_idx:idx) < values(start_idx))
            trough_idx = idx;
            return
        end
    end

    [~, local_idx] = min(values(start_idx:end));
    trough_idx = start_idx + local_idx - 1;
end

function clipped = clip_indices(indices, max_idx)
    clipped = indices(indices >= 1 & indices <= max_idx);
    clipped = unique(clipped, 'stable');
end

function [peak_rate_hz, peak_idx, chosen_window] = first_event_peak(mean_rate, onset_window, offset_window)
    [onset_peak, onset_idx] = local_peak(mean_rate, onset_window);
    [offset_peak, offset_idx] = local_peak(mean_rate, offset_window);

    if onset_peak >= offset_peak
        peak_rate_hz = onset_peak;
        peak_idx = onset_idx;
        chosen_window = onset_window;
    else
        peak_rate_hz = offset_peak;
        peak_idx = offset_idx;
        chosen_window = offset_window;
    end

    if ~isfinite(peak_rate_hz)
        peak_rate_hz = NaN;
        peak_idx = NaN;
        chosen_window = [];
    end
end

function [peak_value, peak_idx] = local_peak(mean_rate, window_indices)
    if isempty(window_indices)
        peak_value = NaN;
        peak_idx = NaN;
        return
    end

    [peak_value, local_idx] = max(mean_rate(window_indices));
    peak_idx = window_indices(local_idx);
end

function plot_stimulus_event_windows(stimulus_event_info, output_name)
    f = figure('Position', [50, 50, 1200, 450]);
    hold on

    x_ms = stimulus_event_info.bin_centers_ms;
    y = stimulus_event_info.envelope_bins;
    onset_x = x_ms(stimulus_event_info.onset_bins);
    offset_x = x_ms(stimulus_event_info.offset_bins);

    if ~isempty(onset_x)
        patch([onset_x(1), onset_x(end), onset_x(end), onset_x(1)], [0, 0, 1.05, 1.05], ...
              [0.85, 0.92, 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.35);
    end

    if ~isempty(offset_x)
        patch([offset_x(1), offset_x(end), offset_x(end), offset_x(1)], [0, 0, 1.05, 1.05], ...
              [1.00, 0.90, 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.35);
    end

    plot(x_ms, y, 'k', 'LineWidth', 1.8)
    xline(x_ms(stimulus_event_info.start_idx), '--', 'Start');
    xline(x_ms(stimulus_event_info.first_peak_idx), '--', '1st peak');
    xline(x_ms(stimulus_event_info.first_trough_idx), '--', '1st trough');

    xlabel('Time (ms)')
    ylabel('Normalized stimulus envelope')
    title('Target1 first onset and first offset windows used for peak features')
    ylim([0, 1.05])
    box on
    grid on
    exportgraphics(f, output_name, 'ContentType', 'vector')
end

function mean_r = mean_pairwise_trial_corr(trial_bins)
    n_trials_local = size(trial_bins, 1);
    pair_rs = nan(nchoosek(n_trials_local, 2), 1);
    count = 0;

    for first_idx = 1:(n_trials_local - 1)
        for second_idx = (first_idx + 1):n_trials_local
            pair_r = safe_corr(trial_bins(first_idx, :).', trial_bins(second_idx, :).');
            if isfinite(pair_r)
                count = count + 1;
                pair_rs(count) = pair_r;
            end
        end
    end

    if count == 0
        mean_r = NaN;
    else
        mean_r = mean(pair_rs(1:count));
    end
end

function focus = get_focus_index(tuning_value)
    tuning_text = lower(char(string(tuning_value)));

    if contains(tuning_text, 'contra')
        focus = 1;
    elseif contains(tuning_text, '45')
        focus = 2;
    elseif contains(tuning_text, 'center')
        focus = 3;
    elseif contains(tuning_text, 'ipsi')
        focus = 4;
    else
        focus = 1;
    end
end

function plot_feature_grid(Y, features_all, feature_display_names, axis1_corrs_all, axis2_corrs_all, ...
                           n_cells, zoom_xlim, zoom_ylim, output_name)
    f = figure('Position', [50, 50, 2600, 1200]);
    t = tiledlayout(2, 5, 'Padding', 'compact', 'TileSpacing', 'compact');

    for feature_idx = 1:numel(feature_display_names)
        ax = nexttile(t, feature_idx);
        hold(ax, 'on')
        draw_pair_lines(ax, Y, n_cells, [0.85, 0.85, 0.85]);

        values = features_all(:, feature_idx);
        valid_mask = isfinite(values);
        scatter(ax, Y(valid_mask, 1), Y(valid_mask, 2), 26, values(valid_mask), 'filled');

        xlabel(ax, 'MDS 1')
        ylabel(ax, 'MDS 2')
        title(ax, sprintf('%s\nr_1 = %.2f | r_2 = %.2f', ...
            feature_display_names{feature_idx}, axis1_corrs_all(feature_idx), axis2_corrs_all(feature_idx)))
        grid(ax, 'on')
        box(ax, 'on')
        colorbar(ax)

        if ~isempty(zoom_xlim)
            xlim(ax, zoom_xlim)
        end
        if ~isempty(zoom_ylim)
            ylim(ax, zoom_ylim)
        end
    end

    if isempty(zoom_xlim)
        sgtitle(t, 'Joint MDS colored by raster/PSTH features')
    else
        sgtitle(t, 'Joint MDS colored by raster/PSTH features (zoomed)')
    end

    exportgraphics(f, output_name, 'ContentType', 'vector')
end

function draw_pair_lines(ax, Y, n_cells, line_color)
    for k = 1:n_cells
        plot(ax, [Y(k, 1), Y(k + n_cells, 1)], ...
                 [Y(k, 2), Y(k + n_cells, 2)], ...
                 'Color', line_color, 'LineWidth', 0.5);
    end
end

function r = safe_corr(x, y)
    x = x(:);
    y = y(:);
    valid_mask = isfinite(x) & isfinite(y);
    x = x(valid_mask);
    y = y(valid_mask);

    if numel(x) < 2 || all(abs(x - x(1)) < eps) || all(abs(y - y(1)) < eps)
        r = NaN;
        return
    end

    corr_matrix = corrcoef(x, y);
    r = corr_matrix(1, 2);
end
