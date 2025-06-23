clear all; close all; clc;

% --- Simulation Parameters ---
FS = 1000; NFFT = 512; M_windows = 10;
signal_frequencies_nominal = [82, 84, 86, 88, 90, 92, 94, 96];
noise_only_frequencies_nominal = [41, 51, 61, 71, 101, 111, 121, 131]; % For FP/TN in H1
fp_reference_frequencies_nominal = noise_only_frequencies_nominal; % For H0 thresholding

num_signal_freqs = length(signal_frequencies_nominal);
num_noise_only_freqs = length(noise_only_frequencies_nominal);
num_fp_ref_freqs = length(fp_reference_frequencies_nominal);

target_pfa = 0.05;
N_trials_H0 = 5000; % Reduced for quicker demo, increase for accuracy
N_master_trials_H1 = 500; % Reduced for quicker demo
snr_db_per_bin_range = -25:5:5; % Coarser SNR for quicker demo

L_neighbors_lft = 2; % Number of neighbors on each side for LFT

total_samples = NFFT * M_windows;
t_segment = (0:NFFT-1)' / FS;
hanning_win = hann(NFFT);

% --- Helper: Adjust Frequencies and Get Bins (MATLAB indexing) ---
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


metrics_to_run = {'MSC', 'CSM', 'LFT'};
sigma_n_sq_time = 1.0; % Noise variance in time domain
P_n_bin_expected = (2.0 / NFFT) * sigma_n_sq_time; % Expected noise power per bin (one-sided)

for metric_idx = 1:length(metrics_to_run)
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
        if mod(i, round(N_trials_H0/10)) == 0 && N_trials_H0 >=10; fprintf('  H0 trial %d/%d for %s\n', i, N_trials_H0, current_metric_name); end
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
    results_H1_metric = struct();
    for c_name_idx = 1:num_criteria
        crit_name_fn = matlab.lang.makeValidName(criteria_names{c_name_idx}); % Field name
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
        A_signal = sqrt(2 * P_s_bin_target); % Amplitude for one-sided spectrum power

        % Aggregators for this SNR
        counts_PD_agg = zeros(1, num_criteria);
        counts_FP_H1_agg = zeros(1, num_criteria);

        for master_trial = 1:N_master_trials_H1
            n_t_all_windows_master = sqrt(sigma_n_sq_time) * randn(total_samples, 1);
            
            % --- Evaluate FP/TN at noise_only_frequencies for this master trial's noise ---
            Y_fft_noise_master_all_windows = zeros(NFFT/2 + 1, M_windows);
            for k_win = 1:M_windows
                start_idx = (k_win-1)*NFFT + 1; end_idx = k_win*NFFT;
                window_data_noise = n_t_all_windows_master(start_idx:end_idx);
                fft_c_noise = fft(window_data_noise .* hanning_win, NFFT);
                Y_fft_noise_master_all_windows(:, k_win) = fft_c_noise(1 : NFFT/2 + 1);
            end
            
            metric_spectrum_noise_master = [];
            switch current_metric_name
                case 'MSC'
                    metric_spectrum_noise_master = msc_fft(Y_fft_noise_master_all_windows, M_windows);
                case 'CSM'
                    metric_spectrum_noise_master = csm_fft_local(Y_fft_noise_master_all_windows, M_windows);
                case 'LFT'
                    metric_spectrum_noise_master = lft_fft_local(Y_fft_noise_master_all_windows, M_windows, L_neighbors_lft);
            end

            metric_at_noise_only_bins_master = metric_spectrum_noise_master(noise_only_bin_indices_matlab);
            fp_results_for_noise_bins_master = (metric_at_noise_only_bins_master > threshold_metric);
            
            num_fp_in_noise_master_trial = sum(fp_results_for_noise_bins_master);
            counts_FP_H1_agg(1) = counts_FP_H1_agg(1) + num_fp_in_noise_master_trial; % Media
            if num_fp_in_noise_master_trial > 0; counts_FP_H1_agg(2) = counts_FP_H1_agg(2) + 1; end % Se Alguma
            if num_fp_in_noise_master_trial > (num_noise_only_freqs / 2); counts_FP_H1_agg(3) = counts_FP_H1_agg(3) + 1; end % Se Maioria
            if num_fp_in_noise_master_trial == num_noise_only_freqs; counts_FP_H1_agg(4) = counts_FP_H1_agg(4) + 1; end % Se Todas

            % --- Evaluate PD/FN at signal_frequencies ---
            individual_signal_detections_this_master_trial = zeros(1, num_signal_freqs);
            for sig_f_idx = 1:num_signal_freqs
                current_f_signal_adj = signal_frequencies_adj(sig_f_idx);
                current_signal_actual_bin = signal_bin_indices_matlab(sig_f_idx);
                
                s_t_single_window = A_signal * sin(2 * pi * current_f_signal_adj * t_segment);
                s_t_all_windows = repmat(s_t_single_window, M_windows, 1);
                x_t_all_windows = s_t_all_windows + n_t_all_windows_master; % Signal + same master noise
                
                Y_fft_H1_all_windows = zeros(NFFT/2 + 1, M_windows);
                for k_win = 1:M_windows
                    start_idx = (k_win-1)*NFFT + 1; end_idx = k_win*NFFT;
                    window_data = x_t_all_windows(start_idx:end_idx);
                    fft_c = fft(window_data .* hanning_win, NFFT);
                    Y_fft_H1_all_windows(:, k_win) = fft_c(1 : NFFT/2 + 1);
                end
                
                metric_spectrum_H1 = [];
                switch current_metric_name
                    case 'MSC'
                        metric_spectrum_H1 = msc_fft(Y_fft_H1_all_windows, M_windows);
                    case 'CSM'
                        metric_spectrum_H1 = csm_fft_local(Y_fft_H1_all_windows, M_windows);
                    case 'LFT'
                        metric_spectrum_H1 = lft_fft_local(Y_fft_H1_all_windows, M_windows, L_neighbors_lft);
                end
                
                current_metric_at_signal_bin = NaN;
                if ~isempty(metric_spectrum_H1) && length(metric_spectrum_H1) >= current_signal_actual_bin && current_signal_actual_bin > 0
                     current_metric_at_signal_bin = metric_spectrum_H1(current_signal_actual_bin);
                end
                if ~isnan(current_metric_at_signal_bin) && current_metric_at_signal_bin > threshold_metric
                    individual_signal_detections_this_master_trial(sig_f_idx) = 1;
                end
            end
            num_detected_signals_in_master_trial = sum(individual_signal_detections_this_master_trial);
            counts_PD_agg(1) = counts_PD_agg(1) + num_detected_signals_in_master_trial; % Media
            if num_detected_signals_in_master_trial > 0; counts_PD_agg(2) = counts_PD_agg(2) + 1; end % Se Alguma
            if num_detected_signals_in_master_trial > (num_signal_freqs / 2); counts_PD_agg(3) = counts_PD_agg(3) + 1; end % Se Maioria
            if num_detected_signals_in_master_trial == num_signal_freqs; counts_PD_agg(4) = counts_PD_agg(4) + 1; end % Se Todas
        end

        % Calculate and store rates for this SNR
        n_PD_denominators = [N_master_trials_H1 * num_signal_freqs, N_master_trials_H1, N_master_trials_H1, N_master_trials_H1];
        n_FP_H1_denominators = [N_master_trials_H1 * num_noise_only_freqs, N_master_trials_H1, N_master_trials_H1, N_master_trials_H1];

        for c_idx = 1:num_criteria
            crit_name_fn = matlab.lang.makeValidName(criteria_names{c_idx});
            
            [pd_val, pd_ci] = wilson_score_interval(counts_PD_agg(c_idx), n_PD_denominators(c_idx), 0.95);
            results_H1_metric.(crit_name_fn).PD(snr_idx) = pd_val;
            results_H1_metric.(crit_name_fn).PD_ic_lower(snr_idx) = pd_ci(1);
            results_H1_metric.(crit_name_fn).PD_ic_upper(snr_idx) = pd_ci(2);
            results_H1_metric.(crit_name_fn).FN_Rate(snr_idx) = 1 - pd_val;
            results_H1_metric.(crit_name_fn).FN_ic_lower(snr_idx) = 1 - pd_ci(2);
            results_H1_metric.(crit_name_fn).FN_ic_upper(snr_idx) = 1 - pd_ci(1);

            [fp_h1_val, fp_h1_ci] = wilson_score_interval(counts_FP_H1_agg(c_idx), n_FP_H1_denominators(c_idx), 0.95);
            results_H1_metric.(crit_name_fn).FP_Rate_H1(snr_idx) = fp_h1_val;
            results_H1_metric.(crit_name_fn).FP_H1_ic_lower(snr_idx) = fp_h1_ci(1);
            results_H1_metric.(crit_name_fn).FP_H1_ic_upper(snr_idx) = fp_h1_ci(2);
            results_H1_metric.(crit_name_fn).TN_Rate_H1(snr_idx) = 1 - fp_h1_val;
            results_H1_metric.(crit_name_fn).TN_H1_ic_lower(snr_idx) = 1 - fp_h1_ci(2);
            results_H1_metric.(crit_name_fn).TN_H1_ic_upper(snr_idx) = 1 - fp_h1_ci(1);
        end
        fprintf('  SNR=%.1fdB (%s): PD_M=%.3f, FP_H1_M=%.4f\n', ...
            current_snr_db, current_metric_name, results_H1_metric.Media.PD(snr_idx), results_H1_metric.Media.FP_Rate_H1(snr_idx));
    end
    time_h1 = toc(tic_h1);
    fprintf('Time total H1 for %s: %.2fs\n', current_metric_name, time_h1);

    % --- Part 3: Plotting for current_metric_name ---
    fig_h = figure('Position', [50, 50, 1200, 800], 'Name', sprintf('Metric: %s', current_metric_name));
    plot_styles = {'-o', '-s', '-^', '-d'};
    colors = lines(num_criteria);

    % Subplot 1: PD
    ax1 = subplot(2,2,1); hold(ax1, 'on');
    for c_idx = 1:num_criteria
        crit_name_fn = matlab.lang.makeValidName(criteria_names{c_idx});
        plot(ax1, snr_db_per_bin_range, results_H1_metric.(crit_name_fn).PD, plot_styles{c_idx}, ...
             'Color', colors(c_idx,:), 'LineWidth', 1.5, 'DisplayName', ['PD ' criteria_names{c_idx}]);
    end
    yline(ax1, target_pfa, 'k--', 'DisplayName', sprintf('PFA Alvo = %.2f', target_pfa));
    title(ax1, sprintf('PD vs. SNR (PFA Emp H0: %.3f)', pfa_empirical_global_H0_metric));
    ylabel(ax1, 'PD (TP Rate)'); xlabel(ax1, 'SNR por Bin (dB)'); grid(ax1, 'on'); legend(ax1, 'show', 'Location','SouthEast'); ylim(ax1, [-0.05, 1.05]);

    % Subplot 2: FN Rate
    ax2 = subplot(2,2,2); hold(ax2, 'on');
    for c_idx = 1:num_criteria
        crit_name_fn = matlab.lang.makeValidName(criteria_names{c_idx});
        plot(ax2, snr_db_per_bin_range, results_H1_metric.(crit_name_fn).FN_Rate, plot_styles{c_idx}, ...
             'Color', colors(c_idx,:), 'LineWidth', 1.5, 'DisplayName', ['FN ' criteria_names{c_idx}]);
    end
    title(ax2, 'FN Rate vs. SNR');
    ylabel(ax2, 'FN Rate'); xlabel(ax2, 'SNR por Bin (dB)'); grid(ax2, 'on'); legend(ax2, 'show', 'Location','NorthEast'); ylim(ax2, [-0.05, 1.05]);

    % Subplot 3: FP Rate (H1)
    ax3 = subplot(2,2,3); hold(ax3, 'on');
    for c_idx = 1:num_criteria
        crit_name_fn = matlab.lang.makeValidName(criteria_names{c_idx});
        plot(ax3, snr_db_per_bin_range, results_H1_metric.(crit_name_fn).FP_Rate_H1, plot_styles{c_idx}, ...
             'Color', colors(c_idx,:), 'LineWidth', 1.5, 'DisplayName', ['FP (H1) ' criteria_names{c_idx}]);
    end
    yline(ax3, target_pfa, 'k--', 'DisplayName', sprintf('PFA Alvo = %.2f', target_pfa));
    yline(ax3, pfa_empirical_global_H0_metric, 'r:', 'DisplayName', sprintf('PFA Emp H0 = %.3f', pfa_empirical_global_H0_metric));
    title(ax3, 'FP Rate (em Freqs Ruído H1) vs. SNR do Sinal');
    ylabel(ax3, 'FP Rate (H1)'); xlabel(ax3, 'SNR por Bin (dB) [do Sinal]'); grid(ax3, 'on'); legend(ax3, 'show', 'Location','NorthEast'); ylim(ax3, [-0.05, 0.5]);

    % Subplot 4: TN Rate (H1)
    ax4 = subplot(2,2,4); hold(ax4, 'on');
    for c_idx = 1:num_criteria
        crit_name_fn = matlab.lang.makeValidName(criteria_names{c_idx});
        plot(ax4, snr_db_per_bin_range, results_H1_metric.(crit_name_fn).TN_Rate_H1, plot_styles{c_idx}, ...
             'Color', colors(c_idx,:), 'LineWidth', 1.5, 'DisplayName', ['TN (H1) ' criteria_names{c_idx}]);
    end
    yline(ax4, tn_rate_global_H0_metric, 'b:', 'DisplayName', sprintf('TN Rate H0 = %.3f', tn_rate_global_H0_metric));
    title(ax4, 'TN Rate (em Freqs Ruído H1) vs. SNR do Sinal');
    ylabel(ax4, 'TN Rate (H1)'); xlabel(ax4, 'SNR por Bin (dB) [do Sinal]'); grid(ax4, 'on'); legend(ax4, 'show', 'Location','SouthEast'); ylim(ax4, [0.5, 1.05]);
    
    sgtitle(fig_h, sprintf('Análise Detecção Multifrequência para %s (PFA Alvo=%.0f%%, M_{win}=%d, L_{neigh}=%d)', ...
        current_metric_name, target_pfa*100, M_windows, L_neighbors_lft), 'FontSize', 14);

end % End of loop over metrics

disp('Simulação com múltiplas métricas concluída.');

% --- Ensure msc_fft and wilson_score_interval are in path or defined in script ---
% function msc = msc_fft(Y_input, M_win) ... (as provided before)
% function [p_hat, ci] = wilson_score_interval(k, n, confidence_level) ... (as provided before)

% --- ORD Functions ---
% msc_fft should be in your path or defined here
% function msc = msc_fft(Y_input, M_win)
%     % Y_input is (num_bins x M_win)
%     if isempty(Y_input) || M_win == 0
%         msc = zeros(size(Y_input,1),1); return;
%     end
%     current_M_win = min(M_win, size(Y_input, 2));
%     if current_M_win == 0
%         msc = zeros(size(Y_input,1),1); return;
%     end
%     if current_M_win == 1 % MSC is 1 for a single window if defined, or 0 if ill-defined
%         msc = ones(size(Y_input,1),1); % Or handle as per specific definition for M=1
%         % A more robust M=1 might be to return 0 or NaN if it implies no coherence
%         % For now, let's assume it's valid to get 1, or rely on M_win > 1 for meaningful MSC
%         return;
%     end
%     Y_to_use = Y_input(:, 1:current_M_win);
%     numerator = abs(sum(Y_to_use, 2)).^2;
%     denominator = current_M_win * sum(abs(Y_to_use).^2, 2);
%     msc = zeros(size(numerator));
%     valid_den = denominator ~= 0;
%     msc(valid_den) = numerator(valid_den) ./ denominator(valid_den);
%     msc(isnan(msc)) = 0; msc(isinf(msc)) = 0;
% end

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
% wilson_score_interval.m should be in your path
