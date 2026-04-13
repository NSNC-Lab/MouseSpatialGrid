clearvars
close all

analysis_dir = 'C:\Users\ipboy\Desktop\SingleChannelFigures\JointMDS_Space';
cd(analysis_dir)

load('50epoch_100ms_and_300ms_loss_All_cells.mat')
load('all_units_info_with_polished_criteria_modified_perf.mat', 'all_data');

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
    'Onset Responsiveness'; ...
    'Temporal Sparsity'; ...
    'Trial Reliability'; ...
    'Offset Responsiveness' ...
    };

[target1, stimulus_fs] = audioread(stimulus_path);
stimulus_event_info = get_stimulus_event_info(target1, stimulus_fs, sample_rate_hz, bin_width, n_bins);

PSTHs_sim = zeros(n_cells, n_bins);
PSTHs_data = zeros(n_cells, n_bins);
features_sim = nan(n_cells, numel(feature_display_names));
features_data = nan(n_cells, numel(feature_display_names));

final_binned_sse = squeeze(losses(end, 2, :, :)).';
[~, best_batch_idx] = min(final_binned_sse, [], 1);

for k = 1:n_cells
    sim_raster = squeeze(output(k, best_batch_idx(k), :, :, :));
    sim_raster = squeeze(sim_raster);
    if isvector(sim_raster)
        sim_raster = reshape(sim_raster, n_trials, []);
    elseif size(sim_raster, 1) ~= n_trials && size(sim_raster, 2) == n_trials
        sim_raster = sim_raster.';
    elseif size(sim_raster, 1) ~= n_trials
        sim_raster = reshape(sim_raster, n_trials, []);
    end
    sim_raster = double(sim_raster);

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

    data_trial_bins = raster_to_trial_bins(data_raster, bin_width);
    PSTHs_data(k, :) = sum(data_trial_bins, 1);
    features_data(k, :) = compute_response_features(data_trial_bins, sample_rate_hz, bin_width, stimulus_event_info);
end

PSTHs_all = [PSTHs_data; PSTHs_sim];
comp_matrix = zeros(2 * n_cells, 2 * n_cells);
for m = 1:(2 * n_cells)
    diffs = PSTHs_all - PSTHs_all(m, :);
    comp_matrix(m, :) = sqrt(sum(diffs .^ 2, 2)).';
end

[Y, stress] = mdscale(comp_matrix, 2, 'criterion', 'metricsstress');
disp("SeparateMDS stress: " + string(stress))

equal_limits = compute_equal_limits(Y, 0.05);
[zoom_x_equal, zoom_y_equal] = equalize_box(zoom_xlim, zoom_ylim);

data_idx = 1:n_cells;
model_idx = (n_cells + 1):(2 * n_cells);
cells = [1:n_cells, 1:n_cells];

features_all = [features_data; features_sim];
axis1_corrs_all = nan(numel(feature_display_names), 1);
axis2_corrs_all = nan(numel(feature_display_names), 1);
axis1_corrs_data = nan(numel(feature_display_names), 1);
axis2_corrs_data = nan(numel(feature_display_names), 1);
axis1_corrs_model = nan(numel(feature_display_names), 1);
axis2_corrs_model = nan(numel(feature_display_names), 1);

for f_idx = 1:numel(feature_display_names)
    axis1_corrs_all(f_idx) = safe_corr(Y(:, 1), features_all(:, f_idx));
    axis2_corrs_all(f_idx) = safe_corr(Y(:, 2), features_all(:, f_idx));
    axis1_corrs_data(f_idx) = safe_corr(Y(data_idx, 1), features_data(:, f_idx));
    axis2_corrs_data(f_idx) = safe_corr(Y(data_idx, 2), features_data(:, f_idx));
    axis1_corrs_model(f_idx) = safe_corr(Y(model_idx, 1), features_sim(:, f_idx));
    axis2_corrs_model(f_idx) = safe_corr(Y(model_idx, 2), features_sim(:, f_idx));
end

feature_summary = table( ...
    string(feature_display_names), ...
    axis1_corrs_all, axis2_corrs_all, ...
    axis1_corrs_data, axis2_corrs_data, ...
    axis1_corrs_model, axis2_corrs_model, ...
    'VariableNames', {'Feature', 'Axis1_r_all', 'Axis2_r_all', 'Axis1_r_data', 'Axis2_r_data', 'Axis1_r_model', 'Axis2_r_model'});

disp(feature_summary)
writetable(feature_summary, 'SeparateMDS_feature_axis_correlations.csv')

save('SeparateMDS_feature_summary.mat', ...
    'Y', 'stress', 'feature_display_names', 'features_data', 'features_sim', ...
    'axis1_corrs_all', 'axis2_corrs_all', ...
    'axis1_corrs_data', 'axis2_corrs_data', ...
    'axis1_corrs_model', 'axis2_corrs_model', ...
    'equal_limits', 'zoom_x_equal', 'zoom_y_equal')

plot_uncolored_joint(Y, cells, n_cells, equal_limits, 'SeparateMDS_joint_uncolored.pdf')
plot_uncolored_subset(Y(data_idx, :), cells(data_idx), equal_limits, 'Data only in joint embedding', 'SeparateMDS_data_uncolored.pdf')
plot_uncolored_subset(Y(model_idx, :), cells(model_idx), equal_limits, 'Model only in joint embedding', 'SeparateMDS_model_uncolored.pdf')

plot_uncolored_joint(Y, cells, n_cells, [zoom_x_equal; zoom_y_equal], 'SeparateMDS_joint_uncolored_zoomed.pdf')
plot_uncolored_subset(Y(data_idx, :), cells(data_idx), [zoom_x_equal; zoom_y_equal], ...
    'Data only in joint embedding (zoomed)', 'SeparateMDS_data_uncolored_zoomed.pdf')
plot_uncolored_subset(Y(model_idx, :), cells(model_idx), [zoom_x_equal; zoom_y_equal], ...
    'Model only in joint embedding (zoomed)', 'SeparateMDS_model_uncolored_zoomed.pdf')

plot_feature_grid_joint(Y, features_all, feature_display_names, axis1_corrs_all, axis2_corrs_all, ...
    n_cells, equal_limits, 'SeparateMDS_feature_grid_joint.pdf')
plot_feature_grid_subset(Y(data_idx, :), features_data, feature_display_names, axis1_corrs_data, axis2_corrs_data, ...
    equal_limits, 'Data-only feature views in joint embedding', 'SeparateMDS_feature_grid_data.pdf')
plot_feature_grid_subset(Y(model_idx, :), features_sim, feature_display_names, axis1_corrs_model, axis2_corrs_model, ...
    equal_limits, 'Model-only feature views in joint embedding', 'SeparateMDS_feature_grid_model.pdf')

plot_feature_grid_joint(Y, features_all, feature_display_names, axis1_corrs_all, axis2_corrs_all, ...
    n_cells, [zoom_x_equal; zoom_y_equal], 'SeparateMDS_feature_grid_joint_zoomed.pdf')
plot_feature_grid_subset(Y(data_idx, :), features_data, feature_display_names, axis1_corrs_data, axis2_corrs_data, ...
    [zoom_x_equal; zoom_y_equal], 'Data-only feature views in joint embedding (zoomed)', 'SeparateMDS_feature_grid_data_zoomed.pdf')
plot_feature_grid_subset(Y(model_idx, :), features_sim, feature_display_names, axis1_corrs_model, axis2_corrs_model, ...
    [zoom_x_equal; zoom_y_equal], 'Model-only feature views in joint embedding (zoomed)', 'SeparateMDS_feature_grid_model_zoomed.pdf')

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
    onset_window = clip_indices(stimulus_event_info.first_onset_bins, n_bins_local);
    offset_window = clip_indices(stimulus_event_info.first_offset_bins, n_bins_local);
    [peak_rate_hz, peak_idx, chosen_window] = first_event_peak(mean_rate, onset_window, offset_window);
    onset_responsiveness = window_responsiveness_index(mean_rate, clip_indices(stimulus_event_info.first_onset_bins, n_bins_local));
    offset_responsiveness = window_responsiveness_index(mean_rate, clip_indices(stimulus_event_info.first_offset_bins, n_bins_local));

    if peak_rate_hz <= 0 || isempty(chosen_window)
        peak_latency_ms = NaN;
        onset_latency_ms = NaN;
        half_max_width_ms = NaN;
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
    end

    total_rate = sum(mean_rate);
    if total_rate > 0
        temporal_center_ms = sum(bin_centers_ms .* mean_rate) / total_rate;
    else
        temporal_center_ms = NaN;
    end

    if all(mean_rate == 0)
        temporal_sparsity = NaN;
    else
        numerator = 1 - ((mean(mean_rate) ^ 2) / (mean(mean_rate .^ 2) + eps));
        denominator = 1 - (1 / n_bins_local);
        temporal_sparsity = numerator / (denominator + eps);
    end

    trial_reliability = mean_pairwise_trial_corr(trial_bins);

    feature_vector = [ ...
        mean_rate_hz, peak_rate_hz, peak_latency_ms, onset_latency_ms, half_max_width_ms, ...
        temporal_center_ms, onset_responsiveness, temporal_sparsity, trial_reliability, offset_responsiveness ...
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
    active_mask = envelope_bins > threshold;
    start_idx = find(active_mask, 1, 'first');
    if isempty(start_idx)
        start_idx = 1;
    end
    end_idx = find(active_mask, 1, 'last');
    if isempty(end_idx)
        end_idx = n_bins;
    end

    first_peak_idx = find_first_local_max(envelope_bins, start_idx);
    if first_peak_idx <= start_idx
        first_peak_idx = min(n_bins, start_idx + 1);
    end

    first_trough_idx = find_first_local_min(envelope_bins, first_peak_idx);
    if first_trough_idx <= first_peak_idx
        first_trough_idx = min(n_bins, first_peak_idx + 1);
    end

    [all_onset_bins, all_offset_bins] = build_rise_fall_masks(envelope_bins, start_idx, end_idx);

    bin_ms = 1000 * (bin_width / response_fs);
    stimulus_event_info = struct();
    stimulus_event_info.envelope_bins = envelope_bins;
    stimulus_event_info.bin_centers_ms = ((1:n_bins) - 0.5) * bin_ms;
    stimulus_event_info.start_idx = start_idx;
    stimulus_event_info.end_idx = end_idx;
    stimulus_event_info.first_peak_idx = first_peak_idx;
    stimulus_event_info.first_trough_idx = first_trough_idx;
    stimulus_event_info.first_onset_bins = start_idx:first_peak_idx;
    stimulus_event_info.first_offset_bins = first_peak_idx:first_trough_idx;
    stimulus_event_info.onset_bins = all_onset_bins;
    stimulus_event_info.offset_bins = all_offset_bins;
end

function peak_idx = find_first_local_max(values, start_idx)
    peak_idx = start_idx;
    for idx = (start_idx + 1):(numel(values) - 1)
        if values(idx) >= values(idx - 1) && values(idx) >= values(idx + 1) && any(values(start_idx:idx) > values(start_idx))
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
        if values(idx) <= values(idx - 1) && values(idx) <= values(idx + 1) && any(values(start_idx:idx) < values(start_idx))
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

function [onset_bins, offset_bins] = build_rise_fall_masks(envelope_bins, start_idx, end_idx)
    onset_bins = [];
    offset_bins = [];

    if end_idx <= start_idx
        onset_bins = start_idx;
        return
    end

    delta = diff(envelope_bins(start_idx:end_idx));
    if isempty(delta)
        onset_bins = start_idx;
        return
    end

    slope_threshold = max(1e-4, 0.01 * max(abs(delta)));
    local_state = zeros(size(delta));
    local_state(delta > slope_threshold) = 1;
    local_state(delta < -slope_threshold) = -1;
    local_state = fill_local_state(local_state);

    if all(local_state == 0)
        onset_bins = start_idx:end_idx;
        return
    end

    bin_state = zeros(1, end_idx - start_idx + 1);
    bin_state(1:end-1) = local_state;
    bin_state(end) = local_state(end);

    onset_bins = (start_idx - 1) + find(bin_state > 0);
    offset_bins = (start_idx - 1) + find(bin_state < 0);

    onset_bins = unique(onset_bins, 'stable');
    offset_bins = unique(offset_bins, 'stable');
end

function local_state = fill_local_state(local_state)
    first_nonzero = find(local_state ~= 0, 1, 'first');
    if isempty(first_nonzero)
        return
    end
    local_state(1:first_nonzero-1) = local_state(first_nonzero);
    for idx = (first_nonzero + 1):numel(local_state)
        if local_state(idx) == 0
            local_state(idx) = local_state(idx - 1);
        end
    end
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

function responsiveness = window_responsiveness_index(mean_rate, focus_window)
    focus_window = unique(focus_window, 'stable');
    focus_window = focus_window(focus_window >= 1 & focus_window <= numel(mean_rate));

    if isempty(focus_window)
        responsiveness = NaN;
        return
    end

    focus_mask = false(size(mean_rate));
    focus_mask(focus_window) = true;
    background_mask = ~focus_mask;
    if ~any(background_mask)
        responsiveness = NaN;
        return
    end

    focus_rate = mean(mean_rate(focus_mask));
    background_rate = mean(mean_rate(background_mask));
    responsiveness = (focus_rate - background_rate) / (focus_rate + background_rate + eps);
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

function plot_uncolored_joint(Y, cells, n_cells, axis_limits, output_name)
    f = figure('Position', [50, 50, 1000, 1000]);
    hold on
    draw_pair_lines(gca, Y, n_cells, [0.85, 0.85, 0.85]);
    pair_handle = plot(nan, nan, '-', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 0.8);
    data_handle = scatter(Y(1:n_cells, 1), Y(1:n_cells, 2), 28, [0.10, 0.35, 0.85], 'o', 'filled');
    model_handle = scatter(Y((n_cells + 1):end, 1), Y((n_cells + 1):end, 2), 34, [0.90, 0.35, 0.15], 's', 'filled');
    xlabel('MDS 1')
    ylabel('MDS 2')
    title('Joint embedding with matched square limits')
    legend([pair_handle, data_handle, model_handle], {'Pairs', 'Data', 'Model'}, ...
        'Location', 'northoutside', 'Orientation', 'horizontal')
    apply_axis_limits(gca, axis_limits)
    labels = string(cells);
    text(Y(:, 1), Y(:, 2), labels, 'FontSize', 8, 'Clipping', 'on')
    exportgraphics(f, output_name, 'ContentType', 'vector')
end

function plot_uncolored_subset(Y_subset, labels_subset, axis_limits, plot_title, output_name)
    f = figure('Position', [50, 50, 1000, 1000]);
    scatter(Y_subset(:, 1), Y_subset(:, 2), 28, 'k', 'filled');
    xlabel('MDS 1')
    ylabel('MDS 2')
    title(plot_title)
    apply_axis_limits(gca, axis_limits)
    labels = string(labels_subset);
    text(Y_subset(:, 1), Y_subset(:, 2), labels, 'FontSize', 8, 'Clipping', 'on')
    grid on
    box on
    exportgraphics(f, output_name, 'ContentType', 'vector')
end

function plot_feature_grid_joint(Y, features_all, feature_display_names, axis1_corrs, axis2_corrs, n_cells, axis_limits, output_name)
    f = figure('Position', [50, 50, 2600, 1200]);
    t = tiledlayout(2, 5, 'Padding', 'compact', 'TileSpacing', 'compact');

    for feature_idx = 1:numel(feature_display_names)
        ax = nexttile(t, feature_idx);
        hold(ax, 'on')
        draw_pair_lines(ax, Y, n_cells, [0.87, 0.87, 0.87]);
        values = features_all(:, feature_idx);
        valid_mask = isfinite(values);

        data_valid = valid_mask(1:n_cells);
        model_valid = valid_mask((n_cells + 1):end);
        data_rows = find(data_valid);
        model_rows = n_cells + find(model_valid);
        scatter(ax, Y(data_rows, 1), Y(data_rows, 2), 24, values(data_rows), 'o', 'filled');
        scatter(ax, Y(model_rows, 1), Y(model_rows, 2), ...
            32, values(model_rows), 's', 'filled', 'MarkerEdgeColor', [0.15, 0.15, 0.15], 'LineWidth', 0.35);

        xlabel(ax, 'MDS 1')
        ylabel(ax, 'MDS 2')
        title(ax, sprintf('%s\nr_1 = %.2f | r_2 = %.2f', ...
            feature_display_names{feature_idx}, axis1_corrs(feature_idx), axis2_corrs(feature_idx)))
        apply_axis_limits(ax, axis_limits)
        colorbar(ax)
    end

    sgtitle(t, 'SeparateMDS joint feature views (circles = data, squares = model)')
    exportgraphics(f, output_name, 'ContentType', 'vector')
end

function plot_feature_grid_subset(Y_subset, feature_matrix, feature_display_names, axis1_corrs, axis2_corrs, axis_limits, plot_title, output_name)
    f = figure('Position', [50, 50, 2600, 1200]);
    t = tiledlayout(2, 5, 'Padding', 'compact', 'TileSpacing', 'compact');

    for feature_idx = 1:numel(feature_display_names)
        ax = nexttile(t, feature_idx);
        hold(ax, 'on')
        values = feature_matrix(:, feature_idx);
        valid_mask = isfinite(values);
        scatter(ax, Y_subset(valid_mask, 1), Y_subset(valid_mask, 2), 28, values(valid_mask), 'filled');
        xlabel(ax, 'MDS 1')
        ylabel(ax, 'MDS 2')
        title(ax, sprintf('%s\nr_1 = %.2f | r_2 = %.2f', ...
            feature_display_names{feature_idx}, axis1_corrs(feature_idx), axis2_corrs(feature_idx)))
        apply_axis_limits(ax, axis_limits)
        colorbar(ax)
    end

    sgtitle(t, plot_title)
    exportgraphics(f, output_name, 'ContentType', 'vector')
end

function draw_pair_lines(ax, Y, n_cells, line_color)
    for k = 1:n_cells
        plot(ax, [Y(k, 1), Y(k + n_cells, 1)], [Y(k, 2), Y(k + n_cells, 2)], 'Color', line_color, 'LineWidth', 0.5);
    end
end

function limits = compute_equal_limits(Y, padding_frac)
    lower = min(Y(:));
    upper = max(Y(:));
    span = upper - lower;
    if span <= 0
        span = 1;
    end
    pad = padding_frac * span;
    limits = [lower - pad, upper + pad];
end

function [x_equal, y_equal] = equalize_box(x_limits, y_limits)
    x_mid = mean(x_limits);
    y_mid = mean(y_limits);
    half_span = 0.5 * max(diff(x_limits), diff(y_limits));
    x_equal = [x_mid - half_span, x_mid + half_span];
    y_equal = [y_mid - half_span, y_mid + half_span];
end

function apply_axis_limits(ax, axis_limits)
    if size(axis_limits, 1) == 2
        xlim(ax, axis_limits(1, :));
        ylim(ax, axis_limits(2, :));
    else
        xlim(ax, axis_limits);
        ylim(ax, axis_limits);
    end
    axis(ax, 'square')
    grid(ax, 'on')
    box(ax, 'on')
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
