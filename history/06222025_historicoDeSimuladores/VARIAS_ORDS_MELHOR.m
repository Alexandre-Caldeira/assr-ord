clear all; close all; clc;

% --- Simulation Parameters ---
FS = 1000; NFFT = 512; M_windows = 10;
signal_frequencies_nominal = [82, 84, 86, 88, 90, 92, 94, 96];
noise_only_frequencies_nominal = [41, 51, 61, 71, 101, 111, 121, 131];
fp_reference_frequencies_nominal = noise_only_frequencies_nominal;

num_signal_freqs = length(signal_frequencies_nominal);
num_noise_only_freqs = length(noise_only_frequencies_nominal);
num_fp_ref_freqs = length(fp_reference_frequencies_nominal);

target_pfa = 0.05;
N_trials_H0 = 5000; % Reduced for demo
N_master_trials_H1 = 5000; % Reduced for demo
snr_db_per_bin_range = -25:5:5; % Coarser for demo

L_neighbors_lft = 2;

total_samples = NFFT * M_windows;
t_segment = (0:NFFT-1)' / FS;
hanning_win = hann(NFFT);

adjust_freqs_get_bins_matlab = @(nominal_freqs, NFFT_val, FS_val) ...
    arrayfun(@(f) deal(round(f * NFFT_val / FS_val) * FS_val / NFFT_val, ...
                       round(f * NFFT_val / FS_val) + 1), nominal_freqs, 'UniformOutput', false);

[signal_frequencies_adj_cell, signal_bin_indices_cell] = adjust_freqs_get_bins_matlab(signal_frequencies_nominal, NFFT, FS);
signal_frequencies_adj = cell2mat(signal_frequencies_adj_cell);
signal_bin_indices_matlab = cell2mat(signal_bin_indices_cell);

[noise_only_frequencies_adj_cell, noise_only_bin_indices_cell] = adjust_freqs_get_bins_matlab(noise_only_frequencies_nominal, NFFT, FS);
noise_only_frequencies_adj = cell2mat(noise_only_frequencies_adj_cell);
noise_only_bin_indices_matlab = cell2mat(noise_only_bin_indices_cell);

[fp_reference_frequencies_adj_cell, fp_ref_bin_indices_cell] = adjust_freqs_get_bins_matlab(fp_reference_frequencies_nominal, NFFT, FS);
fp_reference_frequencies_adj = cell2mat(fp_reference_frequencies_adj_cell);
fp_ref_bin_indices_matlab = cell2mat(fp_ref_bin_indices_cell);

fprintf('Signal Freqs Adjusted: %s\n', mat2str(signal_frequencies_adj,3));
fprintf('Noise (H1) Freqs Adjusted: %s\n', mat2str(noise_only_frequencies_adj,3));
fprintf('FP Ref (H0) Freqs Adjusted: %s\n', mat2str(fp_reference_frequencies_adj,3));

% --- ORD Functions (msc_fft, csm_fft_local, lft_fft_local as previously defined) ---
% (Make sure these functions are defined here or are in your MATLAB path)
% function msc = msc_fft(Y_input, M_win) ...
% function csm_spectrum = csm_fft_local(Y_fft_all_windows, ~) ...
% function lft_spectrum = lft_fft_local(Y_fft_all_windows, ~, L_neigh) ... (use the corrected version from previous response)

% (Assuming wilson_score_interval.m is in path)

metrics_to_run = {'MSC', 'CSM', 'LFT'};
num_metrics = length(metrics_to_run);
aggregated_results_se_alguma = cell(num_metrics, 1); % To store results for consolidated plot

sigma_n_sq_time = 1.0;
P_n_bin_expected = (2.0 / NFFT) * sigma_n_sq_time;

for metric_idx = 1:num_metrics
    current_metric_name = metrics_to_run{metric_idx};
    fprintf('\n--- Processing Metric: %s ---\n', current_metric_name);

    % --- Part 1: H0 - Threshold Determination ---
    fprintf('Calculating threshold for %s (H0)...\n', current_metric_name);
    tic_h0 = tic;
    max_metric_values_H0 = zeros(N_trials_H0, 1);

    for i = 1:N_trials_H0
        noise_all_windows_H0 = sqrt(sigma_n_sq_time) * randn(total_samples, 1);
        Y_fft_H0_all_windows = zeros(NFFT/2 + 1, M_windows);
        for k_win = 1:M_windows
            start_idx = (k_win-1)*NFFT + 1; end_idx = k_win*NFFT;
            window_data = noise_all_windows_H0(start_idx:end_idx);
            fft_coeffs = fft(window_data .* hanning_win, NFFT);
            Y_fft_H0_all_windows(:, k_win) = fft_coeffs(1 : NFFT/2 + 1);
        end

        metric_spectrum_H0 = [];
        switch current_metric_name
            case 'MSC'
                metric_spectrum_H0 = msc_fft(Y_fft_H0_all_windows, M_windows);
            case 'CSM'
                metric_spectrum_H0 = csm_fft_local(Y_fft_H0_all_windows, M_windows);
            case 'LFT'
                metric_spectrum_H0 = lft_fft_local(Y_fft_H0_all_windows, M_windows, L_neighbors_lft);
        end
        
        valid_indices_h0 = fp_ref_bin_indices_matlab(fp_ref_bin_indices_matlab > 0 & fp_ref_bin_indices_matlab <= size(metric_spectrum_H0,1));
        if ~isempty(metric_spectrum_H0) && length(valid_indices_h0) == num_fp_ref_freqs
            metric_at_fp_ref_bins = metric_spectrum_H0(valid_indices_h0);
            max_metric_values_H0(i) = max(metric_at_fp_ref_bins);
        else
            max_metric_values_H0(i) = NaN;
        end
    end
    
    valid_metric_H0 = max_metric_values_H0(~isnan(max_metric_values_H0));
    if isempty(valid_metric_H0); error('Could not calculate H0 metric values for %s.', current_metric_name); end
    
    threshold_metric = prctile(valid_metric_H0, (1 - target_pfa) * 100);
    num_fp_global_H0 = sum(valid_metric_H0 > threshold_metric);
    pfa_empirical_global_H0_metric = num_fp_global_H0 / length(valid_metric_H0);
    tn_rate_global_H0_metric = 1 - pfa_empirical_global_H0_metric;
    time_h0 = toc(tic_h0);

    fprintf('Threshold %s: %.4f (Time H0: %.2fs)\n', current_metric_name, threshold_metric, time_h0);
    fprintf('PFA Emp. Global (H0) for %s: %.4f\n', current_metric_name, pfa_empirical_global_H0_metric);
    fprintf('TN Rate Global (H0) for %s: %.4f\n', current_metric_name, tn_rate_global_H0_metric);

    % --- Part 2: H1 - PD/FN (Signal) & FP/TN (Noise in H1) ---
    criteria_names = {'Media', 'Se Alguma', 'Se Maioria', 'Se Todas'};
    num_criteria = length(criteria_names);
    results_H1_metric = struct(); % Temporary for current metric
    for c_name_idx = 1:num_criteria
        crit_name_fn = matlab.lang.makeValidName(criteria_names{c_name_idx});
        results_H1_metric.(crit_name_fn).PD = zeros(length(snr_db_per_bin_range), 1);
        results_H1_metric.(crit_name_fn).PD_ic_lower = zeros(length(snr_db_per_bin_range), 1);
        results_H1_metric.(crit_name_fn).PD_ic_upper = zeros(length(snr_db_per_bin_range), 1);
        results_H1_metric.(crit_name_fn).FN_Rate = zeros(length(snr_db_per_bin_range), 1);
        results_H1_metric.(crit_name_fn).FN_ic_lower = zeros(length(snr_db_per_bin_range), 1);
        results_H1_metric.(crit_name_fn).FN_ic_upper = zeros(length(snr_db_per_bin_range), 1);
        results_H1_metric.(crit_name_fn).FP_Rate_H1 = zeros(length(snr_db_per_bin_range), 1);
        results_H1_metric.(crit_name_fn).FP_H1_ic_lower = zeros(length(snr_db_per_bin_range), 1);
        results_H1_metric.(crit_name_fn).FP_H1_ic_upper = zeros(length(snr_db_per_bin_range), 1);
        results_H1_metric.(crit_name_fn).TN_Rate_H1 = zeros(length(snr_db_per_bin_range), 1);
        results_H1_metric.(crit_name_fn).TN_H1_ic_lower = zeros(length(snr_db_per_bin_range), 1);
        results_H1_metric.(crit_name_fn).TN_H1_ic_upper = zeros(length(snr_db_per_bin_range), 1);
    end

    fprintf('Calculating PD/FN & FP_H1/TN_H1 for %s...\n', current_metric_name);
    tic_h1 = tic;
    for snr_idx = 1:length(snr_db_per_bin_range)
        current_snr_db = snr_db_per_bin_range(snr_idx);
        current_snr_linear = 10^(current_snr_db / 10);
        P_s_bin_target = current_snr_linear * P_n_bin_expected;
        A_signal = sqrt(2 * P_s_bin_target);

        counts_PD_agg = zeros(1, num_criteria);
        counts_FP_H1_agg = zeros(1, num_criteria);

        for master_trial = 1:N_master_trials_H1
            n_t_all_windows_master = sqrt(sigma_n_sq_time) * randn(total_samples, 1);
            
            Y_fft_noise_master_all_windows = zeros(NFFT/2 + 1, M_windows);
            for k_win = 1:M_windows
                start_idx = (k_win-1)*NFFT + 1; end_idx = k_win*NFFT;
                fft_c_noise = fft(n_t_all_windows_master(start_idx:end_idx) .* hanning_win, NFFT);
                Y_fft_noise_master_all_windows(:, k_win) = fft_c_noise(1 : NFFT/2 + 1);
            end
            
            metric_spectrum_noise_master = [];
            switch current_metric_name
                case 'MSC'; metric_spectrum_noise_master = msc_fft(Y_fft_noise_master_all_windows, M_windows);
                case 'CSM'; metric_spectrum_noise_master = csm_fft_local(Y_fft_noise_master_all_windows, M_windows);
                case 'LFT'; metric_spectrum_noise_master = lft_fft_local(Y_fft_noise_master_all_windows, M_windows, L_neighbors_lft);
            end

            metric_at_noise_only_bins_master = metric_spectrum_noise_master(noise_only_bin_indices_matlab);
            fp_results_for_noise_bins_master = (metric_at_noise_only_bins_master > threshold_metric);
            
            num_fp_in_noise_master_trial = sum(fp_results_for_noise_bins_master);
            counts_FP_H1_agg(1) = counts_FP_H1_agg(1) + num_fp_in_noise_master_trial;
            if num_fp_in_noise_master_trial > 0; counts_FP_H1_agg(2) = counts_FP_H1_agg(2) + 1; end
            if num_fp_in_noise_master_trial > (num_noise_only_freqs / 2); counts_FP_H1_agg(3) = counts_FP_H1_agg(3) + 1; end
            if num_fp_in_noise_master_trial == num_noise_only_freqs; counts_FP_H1_agg(4) = counts_FP_H1_agg(4) + 1; end

            individual_signal_detections_this_master_trial = zeros(1, num_signal_freqs);
            for sig_f_idx = 1:num_signal_freqs
                current_f_signal_adj = signal_frequencies_adj(sig_f_idx);
                current_signal_actual_bin = signal_bin_indices_matlab(sig_f_idx);
                s_t_single_window = A_signal * sin(2 * pi * current_f_signal_adj * t_segment);
                s_t_all_windows = repmat(s_t_single_window, M_windows, 1);
                x_t_all_windows = s_t_all_windows + n_t_all_windows_master;
                
                Y_fft_H1_all_windows = zeros(NFFT/2 + 1, M_windows);
                for k_win = 1:M_windows
                    start_idx = (k_win-1)*NFFT + 1; end_idx = k_win*NFFT;
                    fft_c = fft(x_t_all_windows(start_idx:end_idx) .* hanning_win, NFFT);
                    Y_fft_H1_all_windows(:, k_win) = fft_c(1 : NFFT/2 + 1);
                end
                
                metric_spectrum_H1 = [];
                switch current_metric_name
                    case 'MSC'; metric_spectrum_H1 = msc_fft(Y_fft_H1_all_windows, M_windows);
                    case 'CSM'; metric_spectrum_H1 = csm_fft_local(Y_fft_H1_all_windows, M_windows);
                    case 'LFT'; metric_spectrum_H1 = lft_fft_local(Y_fft_H1_all_windows, M_windows, L_neighbors_lft);
                end
                
                current_metric_at_signal_bin = metric_spectrum_H1(current_signal_actual_bin);
                if ~isnan(current_metric_at_signal_bin) && current_metric_at_signal_bin > threshold_metric
                    individual_signal_detections_this_master_trial(sig_f_idx) = 1;
                end
            end
            num_detected_signals_in_master_trial = sum(individual_signal_detections_this_master_trial);
            counts_PD_agg(1) = counts_PD_agg(1) + num_detected_signals_in_master_trial;
            if num_detected_signals_in_master_trial > 0; counts_PD_agg(2) = counts_PD_agg(2) + 1; end
            if num_detected_signals_in_master_trial > (num_signal_freqs / 2); counts_PD_agg(3) = counts_PD_agg(3) + 1; end
            if num_detected_signals_in_master_trial == num_signal_freqs; counts_PD_agg(4) = counts_PD_agg(4) + 1; end
        end

        n_PD_denominators = [N_master_trials_H1 * num_signal_freqs, N_master_trials_H1, N_master_trials_H1, N_master_trials_H1];
        n_FP_H1_denominators = [N_master_trials_H1 * num_noise_only_freqs, N_master_trials_H1, N_master_trials_H1, N_master_trials_H1];

        for c_idx = 1:num_criteria
            crit_name_fn = matlab.lang.makeValidName(criteria_names{c_idx});
            [pd_val, pd_ci] = wilson_score_interval(counts_PD_agg(c_idx), n_PD_denominators(c_idx), 0.95);
            results_H1_metric.(crit_name_fn).PD(snr_idx) = pd_val;
            results_H1_metric.(crit_name_fn).PD_ic_lower(snr_idx) = pd_ci(1); % Store ICs if you plan to plot them
            results_H1_metric.(crit_name_fn).PD_ic_upper(snr_idx) = pd_ci(2);
            results_H1_metric.(crit_name_fn).FN_Rate(snr_idx) = 1 - pd_val;

            [fp_h1_val, ~] = wilson_score_interval(counts_FP_H1_agg(c_idx), n_FP_H1_denominators(c_idx), 0.95);
            results_H1_metric.(crit_name_fn).FP_Rate_H1(snr_idx) = fp_h1_val;
            results_H1_metric.(crit_name_fn).TN_Rate_H1(snr_idx) = 1 - fp_h1_val;
        end
         fprintf('  SNR=%.1fdB (%s): PD_M=%.3f, FP_H1_M=%.4f\n', ...
            current_snr_db, current_metric_name, results_H1_metric.Media.PD(snr_idx), results_H1_metric.Media.FP_Rate_H1(snr_idx));
    end
    time_h1 = toc(tic_h1);
    fprintf('Time total H1 for %s: %.2fs\n', current_metric_name, time_h1);
    
    % Store "Se Alguma" results for the current metric for consolidated plotting
    idx_se_alguma_crit = find(strcmp(criteria_names, 'Se Alguma'));
    crit_name_fn_se_alguma = matlab.lang.makeValidName(criteria_names{idx_se_alguma_crit});
    
    current_metric_agg_data = struct();
    current_metric_agg_data.name = current_metric_name;
    current_metric_agg_data.pfa_empirical_H0 = pfa_empirical_global_H0_metric; % Store H0 PFA
    current_metric_agg_data.tn_rate_H0 = tn_rate_global_H0_metric; % Store H0 TNR
    
    current_metric_agg_data.PD = results_H1_metric.(crit_name_fn_se_alguma).PD;
    current_metric_agg_data.FN_Rate = results_H1_metric.(crit_name_fn_se_alguma).FN_Rate;
    current_metric_agg_data.FP_Rate_H1 = results_H1_metric.(crit_name_fn_se_alguma).FP_Rate_H1;
    current_metric_agg_data.TN_Rate_H1 = results_H1_metric.(crit_name_fn_se_alguma).TN_Rate_H1;
    
    aggregated_results_se_alguma{metric_idx} = current_metric_agg_data;

    % --- Individual Plotting for current_metric_name (Optional - can be commented out) ---
    % fig_h = figure(...); subplot(2,2,1); ... (as before)
    % disp(['Individual plot for ', current_metric_name, ' generated (optional).']);
    
end % End of loop over metrics

% --- Part 4: Consolidated Plotting for "Se Alguma" Criterion ---
disp('--- Generating Consolidated Plot for "Se Alguma" Criterion ---');
fig_consolidated = figure('Position', [50 50 1200 850], 'Name', 'Consolidated Performance ("Se Alguma")');
plot_styles_metrics = {'-o', '-s', '-^', '-d', '-v', '-<'}; % Styles for different metrics
metric_colors = lines(num_metrics);

% Subplot 1: PD (TP Rate) - "Se Alguma"
ax_pd = subplot(2,2,1); hold(ax_pd, 'on');
for m_idx = 1:num_metrics
    data = aggregated_results_se_alguma{m_idx};
    plot(ax_pd, snr_db_per_bin_range, data.PD, plot_styles_metrics{m_idx}, ...
         'Color', metric_colors(m_idx,:), 'LineWidth', 1.5, 'DisplayName', ['PD ' data.name]);
end
yline(ax_pd, target_pfa, 'k--', 'DisplayName', sprintf('PFA Alvo = %.2f', target_pfa));
title(ax_pd, 'PD (TP Rate) vs. SNR - Criterion: Se Alguma');
ylabel(ax_pd, 'PD (TP Rate)'); xlabel(ax_pd, 'SNR por Bin (dB)'); grid(ax_pd, 'on');
legend(ax_pd, 'show', 'Location','SouthEast'); ylim(ax_pd, [-0.05, 1.05]);

% Subplot 2: FN Rate - "Se Alguma"
ax_fn = subplot(2,2,2); hold(ax_fn, 'on');
for m_idx = 1:num_metrics
    data = aggregated_results_se_alguma{m_idx};
    plot(ax_fn, snr_db_per_bin_range, data.FN_Rate, plot_styles_metrics{m_idx}, ...
         'Color', metric_colors(m_idx,:), 'LineWidth', 1.5, 'DisplayName', ['FN ' data.name]);
end
title(ax_fn, 'FN Rate vs. SNR - Criterion: Se Alguma');
ylabel(ax_fn, 'FN Rate'); xlabel(ax_fn, 'SNR por Bin (dB)'); grid(ax_fn, 'on');
legend(ax_fn, 'show', 'Location','NorthEast'); ylim(ax_fn, [-0.05, 1.05]);

% Subplot 3: FP Rate (H1) - "Se Alguma"
ax_fp_h1 = subplot(2,2,3); hold(ax_fp_h1, 'on');
legend_handles_fp = []; 
for m_idx = 1:num_metrics
    data = aggregated_results_se_alguma{m_idx};
    h = plot(ax_fp_h1, snr_db_per_bin_range, data.FP_Rate_H1, plot_styles_metrics{m_idx}, ...
             'Color', metric_colors(m_idx,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('FP (H1) %s (Emp PFA H0: %.3f)', data.name, data.pfa_empirical_H0));
    legend_handles_fp(end+1) = h; %#ok<AGROW>
end
h_target_fp = yline(ax_fp_h1, target_pfa, 'k--', 'DisplayName', sprintf('PFA Alvo = %.2f', target_pfa));
legend_handles_fp(end+1) = h_target_fp; %#ok<AGROW>
title(ax_fp_h1, 'FP Rate (em Freqs Ruído H1) vs. SNR - Criterion: Se Alguma');
ylabel(ax_fp_h1, 'FP Rate (H1)'); xlabel(ax_fp_h1, 'SNR por Bin (dB) [do Sinal]'); grid(ax_fp_h1, 'on');
legend(ax_fp_h1, legend_handles_fp, 'Location','NorthEast'); ylim(ax_fp_h1, [-0.05, 0.5]);

% Subplot 4: TN Rate (H1) - "Se Alguma"
ax_tn_h1 = subplot(2,2,4); hold(ax_tn_h1, 'on');
legend_handles_tn = [];
for m_idx = 1:num_metrics
    data = aggregated_results_se_alguma{m_idx};
    h = plot(ax_tn_h1, snr_db_per_bin_range, data.TN_Rate_H1, plot_styles_metrics{m_idx}, ...
             'Color', metric_colors(m_idx,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('TN (H1) %s (Emp TNR H0: %.3f)', data.name, data.tn_rate_H0));
    legend_handles_tn(end+1) = h; %#ok<AGROW>
end
% Optional reference line: yline(ax_tn_h1, 1-target_pfa, 'k-.', 'DisplayName', sprintf('1 - PFA Alvo = %.2f', 1-target_pfa));
title(ax_tn_h1, 'TN Rate (em Freqs Ruído H1) vs. SNR - Criterion: Se Alguma');
ylabel(ax_tn_h1, 'TN Rate (H1)'); xlabel(ax_tn_h1, 'SNR por Bin (dB) [do Sinal]'); grid(ax_tn_h1, 'on');
legend(ax_tn_h1, legend_handles_tn, 'Location','SouthEast'); ylim(ax_tn_h1, [0.5, 1.05]);

sgtitle(fig_consolidated, sprintf('Consolidated ORD Performance - Criterion: "Se Alguma" (PFA Alvo=%.0f%%, M_{win}=%d)', ...
    target_pfa*100, M_windows), 'FontSize', 14);

disp('Simulação com múltiplas métricas e plot consolidado concluída.');

function csm_spectrum = csm_fft_local(Y_fft_all_windows, ~) % M_win not explicitly used if Y_fft already epoched
    % Y_fft_all_windows is (num_bins x M_win)
    if isempty(Y_fft_all_windows) || size(Y_fft_all_windows,2) == 0
        csm_spectrum = zeros(size(Y_fft_all_windows,1),1); return;
    end
    phases = angle(Y_fft_all_windows);
    mean_cos = mean(cos(phases), 2);
    mean_sin = mean(sin(phases), 2);
    csm_spectrum = mean_cos.^2 + mean_sin.^2;
    csm_spectrum(isnan(csm_spectrum)) = 0;
end

function lft_spectrum = lft_fft_local(Y_fft_all_windows, ~, L_neigh)
    % Y_fft_all_windows is (num_bins x M_win)
    if isempty(Y_fft_all_windows) || size(Y_fft_all_windows,2) == 0
        lft_spectrum = zeros(size(Y_fft_all_windows,1),1); return;
    end
    power_spectrum_per_window = abs(Y_fft_all_windows).^2;
    avg_power_spectrum = mean(power_spectrum_per_window, 2);
    num_bins = size(avg_power_spectrum, 1);
    lft_spectrum = zeros(num_bins, 1);

    for f_idx = 1:num_bins
        signal_power = avg_power_spectrum(f_idx);
        
        lower_bound = max(1, f_idx - L_neigh);
        upper_bound = min(num_bins, f_idx + L_neigh);
        
        neighbor_indices = [ (lower_bound : f_idx-1), (f_idx+1 : upper_bound) ];
        % Ensure neighbor_indices are valid and unique
        neighbor_indices = unique(neighbor_indices(neighbor_indices >= 1 & neighbor_indices <= num_bins & neighbor_indices ~= f_idx));

        if isempty(neighbor_indices) || length(neighbor_indices) < L_neigh % Need sufficient distinct neighbors
            mean_noise_power = 1e-12; % Avoid division by zero, effectively high LFT
        else
            noise_power_bins = avg_power_spectrum(neighbor_indices);
            mean_noise_power = mean(noise_power_bins);
        end

        if mean_noise_power < 1e-12 % Denominator is too small
            lft_spectrum(f_idx) = signal_power / 1e-12; 
        else
            lft_spectrum(f_idx) = signal_power / mean_noise_power;
        end
    end
    
    % Handle NaNs first
    lft_spectrum(isnan(lft_spectrum)) = 0; 

    % Handle Infs - This is the corrected section
    inf_locations = isinf(lft_spectrum); % Find where Infs are
    if any(inf_locations) % Only proceed if there are Infs
        finite_non_nan_vals = lft_spectrum(~inf_locations & ~isnan(lft_spectrum)); % Values that are NOT Inf and NOT NaN
        
        if isempty(finite_non_nan_vals)
            cap_value_for_inf = 1e6; % Default cap if no other finite, non-NaN values exist
        else
            % Cap at current finite max or 1e6, whichever is larger. Ensure scalar.
            cap_value_for_inf = max(max(finite_non_nan_vals), 1e6); 
        end
        lft_spectrum(inf_locations) = cap_value_for_inf; % Apply the scalar cap value
    end
    
    % The subsequent 'if isempty(lft_spectrum(isinf(lft_spectrum))) && any(isinf(lft_spectrum))'
    % block was not effective and can be removed as this new logic handles Inf capping.
end