clearvars
close all

analysis_dir = 'C:\Users\ipboy\Desktop\SingleChannelFigures\JointMDS_Space';
cd(analysis_dir)

%load('50epoch_100ms_and_300ms_loss_All_cells.mat')
load('LongRunResults.mat')


load('all_units_info_with_polished_criteria_modified_perf.mat', 'all_data');

bin_width = 200;
sample_rate_hz = 10000;
time_len = 29801;
n_cells = size(output, 1);
n_trials = size(output, 3);
n_bins = floor(time_len / bin_width);
zoom_xlim = [-35, -15];
zoom_ylim = [-10, 10];

param_display_names = { ...
    'STRF Gain'; ...
    'STRF Latency'; ...
    'Output Adaptation'; ...
    'On -> ROn gSYN'; ...
    'Off -> ROn gSYN'; ...
    'On -> SOnOff gSYN'; ...
    'Off -> SOnOff gSYN'; ...
    'SOnOff -> ROn gSYN' ...
    };

final_binned_sse = squeeze(losses(end, 2, :, :)).';
[~, best_batch_idx] = min(final_binned_sse, [], 1);

best_params = nan(n_cells, numel(param_display_names));
for k = 1:n_cells
    best_params(k, :) = squeeze(params(end, :, k, best_batch_idx(k)));
end

PSTHs_sim = zeros(n_cells, n_bins);
PSTHs_data = zeros(n_cells, n_bins);

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
    PSTHs_sim(k, :) = sum(raster_to_trial_bins(double(sim_raster), bin_width), 1);

    focus = get_focus_index(all_data(k).tuning_type);
    spike_times = all_data(k).ctrl_tar1_timestamps(:, focus);
    data_raster = zeros(n_trials, time_len);
    for m = 1:n_trials
        stim_mask = (spike_times{m} > 0) & (spike_times{m} < 2.9801);
        trial_indices = round(spike_times{m}(stim_mask) * sample_rate_hz);
        trial_indices = trial_indices(trial_indices >= 0 & trial_indices < time_len);
        data_raster(m, trial_indices + 1) = 1;
    end
    PSTHs_data(k, :) = sum(raster_to_trial_bins(data_raster, bin_width), 1);
end

PSTHs_all = [PSTHs_data; PSTHs_sim];
comp_matrix = zeros(2 * n_cells, 2 * n_cells);
for m = 1:(2 * n_cells)
    diffs = PSTHs_all - PSTHs_all(m, :);
    comp_matrix(m, :) = sqrt(sum(diffs .^ 2, 2)).';
end

[Y, stress] = mdscale(comp_matrix, 2, 'criterion', 'metricsstress');
disp("JointMDSparams stress: " + string(stress))

params_all = [best_params; best_params];
axis1_corrs_all = nan(numel(param_display_names), 1);
axis2_corrs_all = nan(numel(param_display_names), 1);
axis1_corrs_model = nan(numel(param_display_names), 1);
axis2_corrs_model = nan(numel(param_display_names), 1);

for p_idx = 1:numel(param_display_names)
    axis1_corrs_all(p_idx) = safe_corr(Y(:, 1), params_all(:, p_idx));
    axis2_corrs_all(p_idx) = safe_corr(Y(:, 2), params_all(:, p_idx));
    axis1_corrs_model(p_idx) = safe_corr(Y((n_cells + 1):end, 1), best_params(:, p_idx));
    axis2_corrs_model(p_idx) = safe_corr(Y((n_cells + 1):end, 2), best_params(:, p_idx));
end

param_summary = table( ...
    string(param_display_names), ...
    axis1_corrs_all, ...
    axis2_corrs_all, ...
    axis1_corrs_model, ...
    axis2_corrs_model, ...
    'VariableNames', {'Parameter', 'Axis1_r_all', 'Axis2_r_all', 'Axis1_r_model', 'Axis2_r_model'});

disp(param_summary)
writetable(param_summary, 'JointMDS_param_axis_correlations.csv')

cell_index = (1:n_cells).';
best_batch_table = table(cell_index, best_batch_idx(:), ...
    'VariableNames', {'Cell', 'BestBatchIdx'});
for p_idx = 1:numel(param_display_names)
    safe_name = matlab.lang.makeValidName(param_display_names{p_idx});
    best_batch_table.(safe_name) = best_params(:, p_idx);
end
writetable(best_batch_table, 'JointMDS_best_parameter_values.csv')

save('JointMDS_param_summary.mat', ...
    'Y', 'stress', 'param_display_names', 'best_batch_idx', 'best_params', ...
    'axis1_corrs_all', 'axis2_corrs_all', 'axis1_corrs_model', 'axis2_corrs_model')

plot_parameter_grid(Y, best_params, param_display_names, axis1_corrs_all, axis2_corrs_all, ...
    axis1_corrs_model, axis2_corrs_model, n_cells, [], [], 'JointMDS_param_grid.pdf')

plot_parameter_grid(Y, best_params, param_display_names, axis1_corrs_all, axis2_corrs_all, ...
    axis1_corrs_model, axis2_corrs_model, n_cells, zoom_xlim, zoom_ylim, 'JointMDS_param_grid_zoomed.pdf')

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

function plot_parameter_grid(Y, best_params, param_display_names, axis1_corrs_all, axis2_corrs_all, ...
        axis1_corrs_model, axis2_corrs_model, n_cells, zoom_xlim, zoom_ylim, output_name)
    f = figure('Position', [50, 50, 2600, 1100]);
    t = tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

    for p_idx = 1:numel(param_display_names)
        ax = nexttile(t, p_idx);
        hold(ax, 'on')
        draw_pair_lines(ax, Y, n_cells, [0.85, 0.85, 0.85]);

        param_values = best_params(:, p_idx);
        valid_mask = isfinite(param_values);
        scatter(ax, Y(valid_mask, 1), Y(valid_mask, 2), 26, param_values(valid_mask), 'o', 'filled', ...
            'MarkerFaceAlpha', 0.90, 'MarkerEdgeColor', 'none');
        scatter(ax, Y(n_cells + find(valid_mask), 1), Y(n_cells + find(valid_mask), 2), ...
            34, param_values(valid_mask), 's', 'filled', 'MarkerFaceAlpha', 0.90, ...
            'MarkerEdgeColor', [0.15, 0.15, 0.15], 'LineWidth', 0.35);

        xlabel(ax, 'MDS 1')
        ylabel(ax, 'MDS 2')
        title(ax, sprintf('%s\nall: r_1 = %.2f | r_2 = %.2f\nmodel: r_1 = %.2f | r_2 = %.2f', ...
            param_display_names{p_idx}, axis1_corrs_all(p_idx), axis2_corrs_all(p_idx), ...
            axis1_corrs_model(p_idx), axis2_corrs_model(p_idx)))
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
        sgtitle(t, 'Joint MDS colored by best-fit parameter values (circles = data, squares = model)')
    else
        sgtitle(t, 'Joint MDS colored by best-fit parameter values (zoomed; circles = data, squares = model)')
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
