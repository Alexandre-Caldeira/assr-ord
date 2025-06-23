clear all; close all; clc;

% --- Parâmetros da Simulação ---
FS = 1000; NFFT_msc = 512; M_windows = 10;
signal_frequencies = [82, 84, 86, 88, 90, 92, 94, 96];
fp_frequencies = [82, 84, 86, 88, 90, 92, 94, 96]; % Para definir o limiar principal
noise_only_frequencies = [30, 50, 70, 120, 150, 180, 200, 220]; % Para análise secundária de FP/TN

num_signal_freqs = length(signal_frequencies);
num_fp_freqs = length(fp_frequencies);
num_noise_only_freqs = length(noise_only_frequencies);

target_pfa = 0.05;
N_trials_H0 = 10000; 
N_master_trials_H1 = 1000;
snr_db_per_bin_range = -25:2:5;

total_samples_msc = NFFT_msc * M_windows;
t_msc_segment = (0:NFFT_msc-1)' / FS;
hanning_win_msc = hann(NFFT_msc);

% Ajustar frequências e obter bins
fp_bin_indices_matlab = zeros(1, num_fp_freqs);
fp_frequencies_adj = zeros(1, num_fp_freqs);
for i=1:num_fp_freqs
    bin_idx_ideal = round(fp_frequencies(i) * NFFT_msc / FS);
    fp_frequencies_adj(i) = bin_idx_ideal * FS / NFFT_msc;
    fp_bin_indices_matlab(i) = bin_idx_ideal + 1;
end
signal_bin_indices_matlab = zeros(1, num_signal_freqs);
signal_frequencies_adjusted = zeros(1, num_signal_freqs);
for i=1:num_signal_freqs
    bin_idx_ideal = round(signal_frequencies(i) * NFFT_msc / FS);
    signal_frequencies_adjusted(i) = bin_idx_ideal * FS / NFFT_msc;
    signal_bin_indices_matlab(i) = bin_idx_ideal + 1;
end
noise_only_bin_indices_matlab = zeros(1, num_noise_only_freqs);
noise_only_frequencies_adj = zeros(1, num_noise_only_freqs);
for i=1:num_noise_only_freqs
    bin_idx_ideal = round(noise_only_frequencies(i) * NFFT_msc / FS);
    noise_only_frequencies_adj(i) = bin_idx_ideal * FS / NFFT_msc;
    noise_only_bin_indices_matlab(i) = bin_idx_ideal + 1;
end
fprintf('Frequências de FP (limiar) ajustadas: %s\n', mat2str(fp_frequencies_adj,3));
fprintf('Frequências de Sinal ajustadas: %s\n', mat2str(signal_frequencies_adjusted,3));
fprintf('Frequências "Apenas de Ruído" ajustadas: %s\n', mat2str(noise_only_frequencies_adj,3));

% --- Parte 1: Determinação do Limiar e PFA/TN Empírica (H0) ---
fprintf('Calculando limiar da MSC e PFA/TN empírica (H0)...\n');
max_msc_values_H0_at_fp_freqs = zeros(N_trials_H0, 1);
msc_values_at_noise_only_freqs_H0 = zeros(N_trials_H0, num_noise_only_freqs); % Para análise secundária
sigma_n_sq_time = 1.0; 

for i = 1:N_trials_H0
    noise_all_windows_H0 = sqrt(sigma_n_sq_time) * randn(total_samples_msc, 1);
    Y_fft_H0_all_windows = zeros(NFFT_msc/2 + 1, M_windows);
    for k = 1:M_windows
        start_idx = (k-1)*NFFT_msc + 1; end_idx = k*NFFT_msc;
        window_data = noise_all_windows_H0(start_idx:end_idx);
        fft_window = fft(window_data .* hanning_win_msc, NFFT_msc);
        Y_fft_H0_all_windows(:, k) = fft_window(1 : NFFT_msc/2 + 1);
    end
    msc_spectrum_H0 = msc_fft(Y_fft_H0_all_windows, M_windows);
    
    % Para limiar principal
    if ~isempty(msc_spectrum_H0) && all(fp_bin_indices_matlab > 0) && all(fp_bin_indices_matlab <= length(msc_spectrum_H0))
        msc_at_fp_bins = msc_spectrum_H0(fp_bin_indices_matlab);
        max_msc_values_H0_at_fp_freqs(i) = max(msc_at_fp_bins);
    else max_msc_values_H0_at_fp_freqs(i) = NaN; end
    
    % Para análise secundária em "noise_only_frequencies"
    if ~isempty(msc_spectrum_H0) && all(noise_only_bin_indices_matlab > 0) && all(noise_only_bin_indices_matlab <= length(msc_spectrum_H0))
        msc_values_at_noise_only_freqs_H0(i, :) = msc_spectrum_H0(noise_only_bin_indices_matlab);
    else msc_values_at_noise_only_freqs_H0(i, :) = NaN; end
        
    if mod(i, round(N_trials_H0/10)) == 0 && N_trials_H0 >=10; fprintf('  H0 trial %d/%d\n', i, N_trials_H0); end
end

% Limiar e PFA/TN principal
valid_msc_H0_for_threshold = max_msc_values_H0_at_fp_freqs(~isnan(max_msc_values_H0_at_fp_freqs));
if isempty(valid_msc_H0_for_threshold); error('Não foi possível calcular valores de MSC para H0 (limiar).'); end
threshold_msc = prctile(valid_msc_H0_for_threshold, (1 - target_pfa) * 100);
num_fp_H0_principal = sum(valid_msc_H0_for_threshold > threshold_msc);
pfa_empirical_H0_principal = num_fp_H0_principal / length(valid_msc_H0_for_threshold);
tn_rate_H0_principal = 1 - pfa_empirical_H0_principal;
fprintf('Limiar MSC (principal): %.4f\n', threshold_msc);
fprintf('PFA Empírica (principal, sobre fp_frequencies): %.4f\n', pfa_empirical_H0_principal);
fprintf('TN Rate (principal, sobre fp_frequencies): %.4f\n', tn_rate_H0_principal);

% PFA/TN em "noise_only_frequencies"
pfa_per_noise_only_freq = zeros(1, num_noise_only_freqs);
tnr_per_noise_only_freq = zeros(1, num_noise_only_freqs);
for j=1:num_noise_only_freqs
    current_freq_mscs_noise_only = msc_values_at_noise_only_freqs_H0(:,j);
    valid_current_freq_mscs_noise_only = current_freq_mscs_noise_only(~isnan(current_freq_mscs_noise_only));
    if isempty(valid_current_freq_mscs_noise_only); 
        pfa_per_noise_only_freq(j) = NaN; tnr_per_noise_only_freq(j) = NaN;
        continue; 
    end
    fp_count = sum(valid_current_freq_mscs_noise_only > threshold_msc);
    tn_count = sum(valid_current_freq_mscs_noise_only <= threshold_msc);
    pfa_per_noise_only_freq(j) = fp_count / length(valid_current_freq_mscs_noise_only);
    tnr_per_noise_only_freq(j) = tn_count / length(valid_current_freq_mscs_noise_only);
end
avg_pfa_in_noise_regions = nanmean(pfa_per_noise_only_freq);
avg_tnr_in_noise_regions = nanmean(tnr_per_noise_only_freq);
fprintf('--- Análise em Frequências "Apenas de Ruído" (H0) ---\n');
fprintf('PFA Empírica Média em regiões de ruído: %.4f\n', avg_pfa_in_noise_regions);
fprintf('TN Rate Média em regiões de ruído: %.4f\n', avg_tnr_in_noise_regions);

% --- Parte 2: Cálculo da PD e FN Rate (H1) - Mesma lógica da resposta anterior ---
% ... (código para pd_crit e fn_crit permanece o mesmo) ...
pd_crit = struct(); 
fn_crit = struct(); 
criteria_names = {'PD Média', 'PD Se Alguma', 'PD Se Maioria', 'PD Se Todas'};
num_criteria = length(criteria_names);
for c_idx = 1:num_criteria
    pd_crit(c_idx).name = criteria_names{c_idx};
    pd_crit(c_idx).values = zeros(length(snr_db_per_bin_range), 1);
    pd_crit(c_idx).ic_lower = zeros(length(snr_db_per_bin_range), 1);
    pd_crit(c_idx).ic_upper = zeros(length(snr_db_per_bin_range), 1);
    fn_crit(c_idx).name = [criteria_names{c_idx} ' (FN)'];
    fn_crit(c_idx).values = zeros(length(snr_db_per_bin_range), 1);
    fn_crit(c_idx).ic_lower = zeros(length(snr_db_per_bin_range), 1);
    fn_crit(c_idx).ic_upper = zeros(length(snr_db_per_bin_range), 1);
end
P_n_bin_expected = (2/NFFT_msc) * sigma_n_sq_time;
fprintf('Calculando PD e FN Rate para diferentes SNRs e Critérios (H1)...\n');
for snr_idx = 1:length(snr_db_per_bin_range)
    current_snr_db = snr_db_per_bin_range(snr_idx);
    current_snr_linear = 10^(current_snr_db / 10);
    P_s_bin_target = current_snr_linear * P_n_bin_expected;
    A_signal = sqrt(2 * P_s_bin_target);
    detections_count_crit = zeros(1, num_criteria);
    for master_trial = 1:N_master_trials_H1
        individual_detections_this_master_trial = zeros(1, num_signal_freqs);
        for sig_f_idx = 1:num_signal_freqs 
            current_f_signal = signal_frequencies_adjusted(sig_f_idx);
            current_signal_actual_bin = signal_bin_indices_matlab(sig_f_idx);
            s_t_single_window = A_signal * sin(2 * pi * current_f_signal * t_msc_segment);
            s_t_all_windows = repmat(s_t_single_window, M_windows, 1);
            n_t_all_windows = sqrt(sigma_n_sq_time) * randn(total_samples_msc, 1);
            x_t_all_windows = s_t_all_windows + n_t_all_windows;
            Y_fft_H1_all_windows = zeros(NFFT_msc/2 + 1, M_windows);
            for k_win = 1:M_windows
                start_idx = (k_win-1)*NFFT_msc + 1; end_idx = k_win*NFFT_msc;
                window_data = x_t_all_windows(start_idx:end_idx);
                fft_window = fft(window_data .* hanning_win_msc, NFFT_msc);
                Y_fft_H1_all_windows(:, k_win) = fft_window(1 : NFFT_msc/2 + 1);
            end
            msc_spectrum_H1 = msc_fft(Y_fft_H1_all_windows, M_windows);
            current_msc_at_signal_bin = NaN;
            if ~isempty(msc_spectrum_H1) && length(msc_spectrum_H1) >= current_signal_actual_bin && current_signal_actual_bin > 0
                 current_msc_at_signal_bin = msc_spectrum_H1(current_signal_actual_bin);
            end
            if ~isnan(current_msc_at_signal_bin) && current_msc_at_signal_bin > threshold_msc
                individual_detections_this_master_trial(sig_f_idx) = 1;
            end
        end 
        num_detected_in_master_trial = sum(individual_detections_this_master_trial);
        detections_count_crit(1) = detections_count_crit(1) + num_detected_in_master_trial;
        if num_detected_in_master_trial > 0; detections_count_crit(2) = detections_count_crit(2) + 1; end
        if num_detected_in_master_trial > (num_signal_freqs / 2); detections_count_crit(3) = detections_count_crit(3) + 1; end
        if num_detected_in_master_trial == num_signal_freqs; detections_count_crit(4) = detections_count_crit(4) + 1; end
    end 
    total_individual_detection_opportunities = N_master_trials_H1 * num_signal_freqs;
    pd_crit(1).values(snr_idx) = detections_count_crit(1) / total_individual_detection_opportunities;
    [~, pci1_pd] = wilson_score_interval(detections_count_crit(1), total_individual_detection_opportunities, 0.95);
    pd_crit(1).ic_lower(snr_idx) = pci1_pd(1); pd_crit(1).ic_upper(snr_idx) = pci1_pd(2);
    fn_crit(1).values(snr_idx) = 1 - pd_crit(1).values(snr_idx);
    fn_crit(1).ic_lower(snr_idx) = 1 - pd_crit(1).ic_upper(snr_idx); 
    fn_crit(1).ic_upper(snr_idx) = 1 - pd_crit(1).ic_lower(snr_idx);
    for c_idx = 2:num_criteria
        pd_crit(c_idx).values(snr_idx) = detections_count_crit(c_idx) / N_master_trials_H1;
        [~, pci_pd] = wilson_score_interval(detections_count_crit(c_idx), N_master_trials_H1, 0.95);
        pd_crit(c_idx).ic_lower(snr_idx) = pci_pd(1); pd_crit(c_idx).ic_upper(snr_idx) = pci_pd(2);
        fn_crit(c_idx).values(snr_idx) = 1 - pd_crit(c_idx).values(snr_idx);
        fn_crit(c_idx).ic_lower(snr_idx) = 1 - pd_crit(c_idx).ic_upper(snr_idx);
        fn_crit(c_idx).ic_upper(snr_idx) = 1 - pd_crit(c_idx).ic_lower(snr_idx);
    end
    fprintf('  SNR=%.1fdB: PD_M=%.3f, PD_Any=%.3f, PD_Maj=%.3f, PD_All=%.3f\n', ...
            current_snr_db, pd_crit(1).values(snr_idx), pd_crit(2).values(snr_idx), ...
            pd_crit(3).values(snr_idx), pd_crit(4).values(snr_idx));
end


% --- Parte 3: Plotagem Combinada "Matriz de Confiança" vs SNR ---
figure('Position', [100, 100, 900, 700]); % Figura maior
hold on;
plot_styles = {'-o', '-s', '-^', '-d'};
colors = get(gca, 'ColorOrder');
if size(colors,1) < num_criteria; colors = [0 0 1; 1 0 0; 0 1 0; 0.5 0.5 0.5]; end % Fallback colors

% Plot PD (TP Rate) para os 4 critérios
legend_handles_pd = [];
for c_idx = 1:num_criteria
    h = plot(snr_db_per_bin_range, pd_crit(c_idx).values, plot_styles{c_idx}, ...
         'Color', colors(c_idx,:), 'LineWidth', 1.5, ...
         'DisplayName', ['TP Rate (' pd_crit(c_idx).name ')']);
    legend_handles_pd = [legend_handles_pd, h];
    fill([snr_db_per_bin_range, fliplr(snr_db_per_bin_range)], ...
         [pd_crit(c_idx).ic_lower', fliplr(pd_crit(c_idx).ic_upper')], ...
         colors(c_idx,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Plot FN Rate para os 4 critérios (usar mesmas cores, mas estilo de linha diferente)
legend_handles_fn = [];
fn_line_styles = {':o', ':s', ':^', ':d'}; % Diferenciar FN
for c_idx = 1:num_criteria
    h = plot(snr_db_per_bin_range, fn_crit(c_idx).values, fn_line_styles{c_idx}, ...
         'Color', colors(c_idx,:)*0.7, 'LineWidth', 1.5, ... % Cor levemente diferente
         'DisplayName', ['FN Rate (' pd_crit(c_idx).name ')']);
    legend_handles_fn = [legend_handles_fn, h];
    fill([snr_db_per_bin_range, fliplr(snr_db_per_bin_range)], ...
         [fn_crit(c_idx).ic_lower', fliplr(fn_crit(c_idx).ic_upper')], ...
         colors(c_idx,:)*0.7, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Plot Linhas Horizontais para FP e TN Rates (globais e em regiões de ruído)
h_pfa_principal = yline(pfa_empirical_H0_principal, '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, ...
      'DisplayName', sprintf('FP Rate (Global PFA Emp: %.3f)', pfa_empirical_H0_principal));
h_tn_principal = yline(tn_rate_H0_principal, '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, ...
      'DisplayName', sprintf('TN Rate (Global Espec: %.3f)', tn_rate_H0_principal));
      
if ~isnan(avg_pfa_in_noise_regions)
    h_pfa_noise_regions = yline(avg_pfa_in_noise_regions, ':', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, ...
          'DisplayName', sprintf('FP Rate (Média Reg. Ruído: %.3f)', avg_pfa_in_noise_regions));
end
if ~isnan(avg_tnr_in_noise_regions)
    h_tn_noise_regions = yline(avg_tnr_in_noise_regions, ':', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5, ...
          'DisplayName', sprintf('TN Rate (Média Reg. Ruído: %.3f)', avg_tnr_in_noise_regions));
end
      
xlabel('SNR por Bin (dB)');
ylabel('Taxa (Rate)');
title_str = sprintf('Análise Completa de Detecção vs. SNR (M_{win}=%d, %d freqs)', M_windows, num_signal_freqs);
title(title_str);

% Montar legenda
all_legend_handles = [legend_handles_pd, legend_handles_fn, h_pfa_principal, h_tn_principal];
if exist('h_pfa_noise_regions', 'var'); all_legend_handles = [all_legend_handles, h_pfa_noise_regions]; end
if exist('h_tn_noise_regions', 'var'); all_legend_handles = [all_legend_handles, h_tn_noise_regions]; end
legend(all_legend_handles, 'Location', 'EastOutside', 'FontSize', 8);

grid on;
ylim([-0.05, 1.05]);
hold off;

disp('Simulação com plot combinado concluída.');
% --- Funções Auxiliares (msc_fft e wilson_score_interval devem estar no path ou definidas no script) ---