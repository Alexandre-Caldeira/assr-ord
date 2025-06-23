clear all; close all; clc;

% --- Parâmetros da Simulação ---
FS = 1000;                  % Taxa de amostragem (Hz)
NFFT_msc = 512;             % Tamanho da FFT para MSC e para definir SNR por bin
M_windows = 10;             % Número de janelas para MSC
f_signal = 100;             % Frequência do sinal (Hz)

target_pfa = 0.05;          % PFA alvo
N_trials_H0 = 50000;        % Trials para estimar limiar H0
N_trials_H1 = 5000;         % Trials por SNR para estimar PD H1

snr_db_per_bin_range = -30:1:5; % SNR por bin em dB

% Duração total do sinal para M_windows (em amostras)
total_samples_msc = NFFT_msc * M_windows;
t_msc_segment = (0:NFFT_msc-1)' / FS; % Tempo para uma janela da MSC

% Ajustar f_signal para cair exatamente em um bin da NFFT_msc
signal_bin_ideal_idx_msc = round(f_signal * NFFT_msc / FS);
f_signal_adjusted = signal_bin_ideal_idx_msc * FS / NFFT_msc;
fprintf('Frequência do sinal ajustada para: %.2f Hz (bin %d para NFFT=%d)\n', ...
        f_signal_adjusted, signal_bin_ideal_idx_msc, NFFT_msc);
% Índice para arrays MATLAB (1-indexado) para espectro unilateral
% (0 a NFFT/2 -> 1 a NFFT/2+1)
msc_signal_bin_matlab_idx = signal_bin_ideal_idx_msc + 1;


% --- Parte 1: Determinação do Limiar da MSC (H0) ---
fprintf('Calculando limiar da MSC para PFA = %.2f (H0)...\n', target_pfa);
msc_values_H0 = zeros(N_trials_H0, 1);
sigma_n_sq_time = 1.0; % Variância do ruído no tempo fixada em 1

for i = 1:N_trials_H0
    % Gerar ruído para todas as janelas
    noise_all_windows_H0 = sqrt(sigma_n_sq_time) * randn(total_samples_msc, 1);
    
    Y_fft_H0_all_windows = zeros(NFFT_msc/2 + 1, M_windows);
    for k = 1:M_windows
        start_idx = (k-1)*NFFT_msc + 1;
        end_idx = k*NFFT_msc;
        window_data = noise_all_windows_H0(start_idx:end_idx);
        % Aplicar janela de Hann para reduzir vazamento (opcional, mas bom)
        fft_window = fft(window_data .* hann(NFFT_msc), NFFT_msc);
        Y_fft_H0_all_windows(:, k) = fft_window(1 : NFFT_msc/2 + 1);
    end
    msc_spectrum_H0 = msc_fft(Y_fft_H0_all_windows, M_windows);
    if ~isempty(msc_spectrum_H0) && length(msc_spectrum_H0) >= msc_signal_bin_matlab_idx
        msc_values_H0(i) = msc_spectrum_H0(msc_signal_bin_matlab_idx);
    else
        msc_values_H0(i) = NaN; % Caso algo dê errado
    end
    if mod(i, round(N_trials_H0/10)) == 0
        fprintf('  H0 trial %d/%d\n', i, N_trials_H0);
    end
end
msc_values_H0 = msc_values_H0(~isnan(msc_values_H0)); % Remover NaNs se houver
threshold_msc = prctile(msc_values_H0, (1 - target_pfa) * 100);
fprintf('Limiar MSC determinado: %.4f\n', threshold_msc);

% --- Parte 2: Cálculo da PD e Intervalo de Confiança (H1) ---
pd_values = zeros(length(snr_db_per_bin_range), 1);
pd_ic_lower = zeros(length(snr_db_per_bin_range), 1);
pd_ic_upper = zeros(length(snr_db_per_bin_range), 1);

% Potência de ruído esperada por bin para NFFT_msc e sigma_n_sq_time=1
P_n_bin_expected = (2/NFFT_msc) * sigma_n_sq_time;

fprintf('Calculando PD para diferentes SNRs (H1)...\n');
for snr_idx = 1:length(snr_db_per_bin_range)
    current_snr_db = snr_db_per_bin_range(snr_idx);
    current_snr_linear = 10^(current_snr_db / 10);
    
    % Potência de sinal necessária no bin do sinal
    P_s_bin_target = current_snr_linear * P_n_bin_expected;
    % Amplitude A da senoide para atingir P_s_bin_target em uma janela
    A_signal = sqrt(2 * P_s_bin_target);
    
    detections_H1 = 0;
    for i = 1:N_trials_H1
        % Gerar sinal senoidal para uma janela
        s_t_single_window = A_signal * sin(2 * pi * f_signal_adjusted * t_msc_segment);
        % Replicar para todas as janelas (sinal coerente entre janelas)
        s_t_all_windows = repmat(s_t_single_window, M_windows, 1);
        
        % Gerar ruído
        n_t_all_windows = sqrt(sigma_n_sq_time) * randn(total_samples_msc, 1);
        x_t_all_windows = s_t_all_windows + n_t_all_windows;
        
        Y_fft_H1_all_windows = zeros(NFFT_msc/2 + 1, M_windows);
        for k = 1:M_windows
            start_idx = (k-1)*NFFT_msc + 1;
            end_idx = k*NFFT_msc;
            window_data = x_t_all_windows(start_idx:end_idx);
            fft_window = fft(window_data .* hann(NFFT_msc), NFFT_msc);
            Y_fft_H1_all_windows(:, k) = fft_window(1 : NFFT_msc/2 + 1);
        end
        msc_spectrum_H1 = msc_fft(Y_fft_H1_all_windows, M_windows);
        
        current_msc_value = NaN;
        if ~isempty(msc_spectrum_H1) && length(msc_spectrum_H1) >= msc_signal_bin_matlab_idx
             current_msc_value = msc_spectrum_H1(msc_signal_bin_matlab_idx);
        end

        if ~isnan(current_msc_value) && current_msc_value > threshold_msc
            detections_H1 = detections_H1 + 1;
        end
    end
    
    pd_values(snr_idx) = detections_H1 / N_trials_H1;
    
    % Intervalo de Confiança de Wilson Score
    [phat, pci] = wilson_score_interval(detections_H1, N_trials_H1, 0.95);
    % phat deve ser igual a pd_values(snr_idx)
    pd_ic_lower(snr_idx) = pci(1);
    pd_ic_upper(snr_idx) = pci(2);
    
    fprintf('  SNR por bin = %.1f dB: PD = %.3f (IC 95%%: [%.3f, %.3f])\n', ...
            current_snr_db, pd_values(snr_idx), pd_ic_lower(snr_idx), pd_ic_upper(snr_idx));
end

% --- Parte 3: Plotagem ---
figure;
% Linha da PD
plot(snr_db_per_bin_range, pd_values, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'PD Média');
hold on;

% Área sombreada para o intervalo de confiança
fill([snr_db_per_bin_range, fliplr(snr_db_per_bin_range)], ...
     [pd_ic_lower', fliplr(pd_ic_upper')], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
     'DisplayName', 'IC 95% para PD');

% Linha do PFA alvo
yline(target_pfa, 'r--', 'LineWidth', 1, ...
      'DisplayName', sprintf('PFA Alvo = %.2f', target_pfa));
% Linha de 1-PFA (idealmente, PD deveria se aproximar disso em SNR alta se não houver outros erros)
% yline(1-target_pfa, 'k:', 'DisplayName', '1 - PFA Alvo');


xlabel('SNR por Bin (dB)');
ylabel('Probabilidade de Detecção (PD)');
title(sprintf('PD da MSC vs. SNR por Bin (PFA Alvo = %.2f%%, M_{windows}=%d)', target_pfa*100, M_windows));
legend('show', 'Location', 'SouthEast');
grid on;
ylim([-0.05, 1.05]);

% --- Funções Auxiliares ---
% (Coloque msc_fft.m e wilson_score_interval.m no seu path ou no final deste script)
% function [ORD] = msc_fft(Y_input, M_win) ... (Definida anteriormente)
% function [phat, pci] = wilson_score_interval(k, n, confidence_level) ... (Definição abaixo)

disp('Simulação concluída.');