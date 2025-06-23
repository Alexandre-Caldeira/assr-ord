clear all; close all; clc;

% --- Parâmetros da Simulação ---
FS = 1000;                  % Taxa de amostragem (Hz)
NFFT_msc = 512;             % Tamanho da FFT para MSC e para definir SNR por bin
M_windows = 10;             % Número de janelas para MSC

signal_frequencies = [82, 84, 86, 88, 90, 92, 94, 96];
fp_frequencies = [82, 84, 86, 88, 90, 92, 94, 96];
num_signal_freqs = length(signal_frequencies);
num_fp_freqs = length(fp_frequencies);

target_pfa = 0.05;
N_trials_H0 = 10000; % Reduzido para tempo de execução mais rápido no exemplo
N_master_trials_H1 = 1000; % Número de "conjuntos de 8 frequências" testados por SNR

snr_db_per_bin_range = -25:2:5; % SNR por bin em dB (ajustado para ver mais variação)

total_samples_msc = NFFT_msc * M_windows;
t_msc_segment = (0:NFFT_msc-1)' / FS;
hanning_win_msc = hann(NFFT_msc);

% Converter frequências para bins e ajustar
fp_bin_indices_matlab = zeros(1, num_fp_freqs);
fp_frequencies_adj = zeros(1, num_fp_freqs); % Renomeado para evitar sobrescrever
for i=1:num_fp_freqs
    bin_idx_ideal = round(fp_frequencies(i) * NFFT_msc / FS);
    fp_frequencies_adj(i) = bin_idx_ideal * FS / NFFT_msc;
    fp_bin_indices_matlab(i) = bin_idx_ideal + 1;
end
fprintf('Frequências de FP ajustadas:\n'); disp(fp_frequencies_adj');

signal_bin_indices_matlab = zeros(1, num_signal_freqs);
signal_frequencies_adjusted = zeros(1, num_signal_freqs);
for i=1:num_signal_freqs
    bin_idx_ideal = round(signal_frequencies(i) * NFFT_msc / FS);
    signal_frequencies_adjusted(i) = bin_idx_ideal * FS / NFFT_msc;
    signal_bin_indices_matlab(i) = bin_idx_ideal + 1;
end
fprintf('Frequências de Sinal ajustadas:\n'); disp(signal_frequencies_adjusted');

% --- Parte 1: Determinação do Limiar da MSC (H0) ---
fprintf('Calculando limiar da MSC para PFA = %.2f (H0)...\n', target_pfa);
max_msc_values_H0_at_fp_freqs = zeros(N_trials_H0, 1);
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
    
    if ~isempty(msc_spectrum_H0) && all(fp_bin_indices_matlab <= length(msc_spectrum_H0)) && all(fp_bin_indices_matlab > 0)
        msc_at_fp_bins = msc_spectrum_H0(fp_bin_indices_matlab);
        max_msc_values_H0_at_fp_freqs(i) = max(msc_at_fp_bins);
    else
        max_msc_values_H0_at_fp_freqs(i) = NaN;
    end
    if mod(i, round(N_trials_H0/10)) == 0 && N_trials_H0 >=10 ; fprintf('  H0 trial %d/%d\n', i, N_trials_H0); end
end
max_msc_values_H0_at_fp_freqs = max_msc_values_H0_at_fp_freqs(~isnan(max_msc_values_H0_at_fp_freqs));
if isempty(max_msc_values_H0_at_fp_freqs)
    error('Não foi possível calcular valores de MSC para H0. Verifique os índices dos bins.');
end
threshold_msc = prctile(max_msc_values_H0_at_fp_freqs, (1 - target_pfa) * 100);
fprintf('Limiar MSC (baseado no max das freqs de FP) determinado: %.4f\n', threshold_msc);

% --- Parte 2: Cálculo da PD para diferentes Critérios (H1) ---
pd_crit = struct(); % Usar estrutura para armazenar resultados
criteria_names = {'PD Média', 'PD Se Alguma', 'PD Se Maioria', 'PD Se Todas'};
num_criteria = length(criteria_names);

for c_idx = 1:num_criteria
    pd_crit(c_idx).name = criteria_names{c_idx};
    pd_crit(c_idx).values = zeros(length(snr_db_per_bin_range), 1);
    pd_crit(c_idx).ic_lower = zeros(length(snr_db_per_bin_range), 1);
    pd_crit(c_idx).ic_upper = zeros(length(snr_db_per_bin_range), 1);
end

P_n_bin_expected = (2/NFFT_msc) * sigma_n_sq_time;

fprintf('Calculando PD para diferentes SNRs e Critérios (H1)...\n');
for snr_idx = 1:length(snr_db_per_bin_range)
    current_snr_db = snr_db_per_bin_range(snr_idx);
    current_snr_linear = 10^(current_snr_db / 10);
    P_s_bin_target = current_snr_linear * P_n_bin_expected;
    A_signal = sqrt(2 * P_s_bin_target);
    
    detections_count_crit = zeros(1, num_criteria); % Contadores para cada critério
    
    for master_trial = 1:N_master_trials_H1
        individual_detections_this_master_trial = zeros(1, num_signal_freqs); % 0 ou 1
        
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
        
        % Critério 1 (PD Média): Soma todas as detecções individuais
        detections_count_crit(1) = detections_count_crit(1) + num_detected_in_master_trial;
        
        % Critério 2: Alguma detectada
        if num_detected_in_master_trial > 0
            detections_count_crit(2) = detections_count_crit(2) + 1;
        end
        
        % Critério 3: Maioria detectada ( > N/2 )
        if num_detected_in_master_trial > (num_signal_freqs / 2)
             detections_count_crit(3) = detections_count_crit(3) + 1;
        end
        
        % Critério 4: Todas detectadas
        if num_detected_in_master_trial == num_signal_freqs
            detections_count_crit(4) = detections_count_crit(4) + 1;
        end
    end 
    
    % Calcular PDs e ICs
    % Para Critério 1 (PD Média), o 'n' do IC é o total de detecções individuais possíveis
    total_individual_detection_opportunities = N_master_trials_H1 * num_signal_freqs;
    pd_crit(1).values(snr_idx) = detections_count_crit(1) / total_individual_detection_opportunities;
    [~, pci1] = wilson_score_interval(detections_count_crit(1), total_individual_detection_opportunities, 0.95);
    pd_crit(1).ic_lower(snr_idx) = pci1(1); pd_crit(1).ic_upper(snr_idx) = pci1(2);

    % Para Critérios 2, 3, 4, o 'n' do IC é N_master_trials_H1
    for c_idx = 2:num_criteria
        pd_crit(c_idx).values(snr_idx) = detections_count_crit(c_idx) / N_master_trials_H1;
        [~, pci] = wilson_score_interval(detections_count_crit(c_idx), N_master_trials_H1, 0.95);
        pd_crit(c_idx).ic_lower(snr_idx) = pci(1); pd_crit(c_idx).ic_upper(snr_idx) = pci(2);
    end
    
    fprintf('  SNR=%.1fdB: PD_M=%.3f, PD_Any=%.3f, PD_Maj=%.3f, PD_All=%.3f\n', ...
            current_snr_db, pd_crit(1).values(snr_idx), pd_crit(2).values(snr_idx), ...
            pd_crit(3).values(snr_idx), pd_crit(4).values(snr_idx));
end

% --- Parte 3: Plotagem ---
figure;
hold on;
plot_styles = {'-o', '-s', '-^', '-d'};
colors = get(gca, 'ColorOrder');
if size(colors,1) < num_criteria
    colors = [0 0 1; 1 0 0; 0 1 0; 0.5 0.5 0.5; 0.75 0 0.75]; % Fallback
end

for c_idx = 1:num_criteria
    plot(snr_db_per_bin_range, pd_crit(c_idx).values, plot_styles{c_idx}, ...
         'Color', colors(c_idx,:), 'LineWidth', 1.5, 'DisplayName', pd_crit(c_idx).name);
    fill([snr_db_per_bin_range, fliplr(snr_db_per_bin_range)], ...
         [pd_crit(c_idx).ic_lower', fliplr(pd_crit(c_idx).ic_upper')], ...
         colors(c_idx,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

yline(target_pfa, 'k--', 'LineWidth', 1, ...
      'DisplayName', sprintf('PFA Alvo = %.2f', target_pfa));
      
xlabel('SNR por Bin (dB)');
ylabel('Probabilidade de Detecção (PD)');
title(sprintf('Comparação de Critérios de Detecção PD vs. SNR (PFA Alvo = %.2f%%, M_{win}=%d, %d freqs)', target_pfa*100, M_windows, num_signal_freqs));
legend('show', 'Location', 'SouthEast');
grid on;
ylim([-0.05, 1.05]);
hold off;

disp('Simulação multifrequência com múltiplos critérios concluída.');

% --- Funções Auxiliares (msc_fft.m e wilson_score_interval.m devem estar no path ou definidas aqui) ---
% (Reutilizar as mesmas funções msc_fft e wilson_score_interval da resposta anterior)