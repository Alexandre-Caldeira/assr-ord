clear all;
close all;
clc;

% --- Parâmetros da Simulação ---
FS = 1000;             % Taxa de amostragem (Hz)
T_window = 0.5;        % Duração de cada janela para FFT/MSC (s)
NFFT = round(T_window * FS);  % Número de pontos da FFT (tamanho da janela)
M_windows = 10;        % Número de janelas para cálculo da MSC (Reduzido para PD ter mais variação com SNR)
signal_freq = 82;      % Frequência do sinal de interesse (Hz)

SNR_dB_range = -50:2:5; % Faixa de SNRs para testar (dB)
SNR_dB_for_spectrum_plot = -10; % dB, para plotar o espectro

% Parâmetros para PD e PFA
target_PFA = 0.05;     % PFA alvo para definir o limiar
N_trials_H0 = 10000;   % Número de realizações para estimar PFA empírica
N_trials_H1 = 2000;    % Número de realizações por SNR para estimar PD

% --- Calcular Limiar Teórico para MSC ---
% Sob H0, MSC ~ Beta(1, M_windows - 1)
threshold_msc = betainv(1 - target_PFA, 1, M_windows - 1);
fprintf('Limiar Teórico da MSC para PFA = %.3f: %.4f (M_windows=%d)\n', target_PFA, threshold_msc, M_windows);

% Duração total do sinal para M_windows
total_duration = M_windows * T_window;
total_samples = round(total_duration * FS);
t = (0:total_samples-1)'/FS; % Vetor de tempo para o sinal completo

% Encontrar o bin de frequência mais próximo de signal_freq
f_axis = (0:NFFT/2)*FS/NFFT;
[~, sig_bin_idx] = min(abs(f_axis - signal_freq));
fprintf('Frequência do sinal: %.1f Hz, Bin da FFT correspondente: %d (%.2f Hz)\n', signal_freq, sig_bin_idx, f_axis(sig_bin_idx));

% --- Estimar PFA Empírica (H0: Somente Ruído) ---
false_alarms_method1 = 0;
false_alarms_method2 = 0;

fprintf('Estimando PFA empírica com %d trials...\n', N_trials_H0);
for trial_idx = 1:N_trials_H0
    % Método 1 (H0): Ruído Gaussiano simples
    noise_H0_m1 = randn(total_samples, 1);
    % Opcional: escalar para variância unitária, embora MSC seja teoricamente invariante à escala global
    % noise_H0_m1 = noise_H0_m1 / std(noise_H0_m1);

    Y_H0_m1_all_windows = zeros(NFFT/2+1, M_windows);
    for k = 1:M_windows
        start_idx = (k-1)*NFFT + 1;
        end_idx = k*NFFT;
        if end_idx > total_samples; continue; end % Proteção
        window_data = noise_H0_m1(start_idx:end_idx);
        fft_window = fft(window_data.*hann(NFFT), NFFT);
        Y_H0_m1_all_windows(:, k) = fft_window(1:NFFT/2+1);
    end
    msc_H0_m1 = msc_fft(Y_H0_m1_all_windows, M_windows);
    if msc_H0_m1(sig_bin_idx) > threshold_msc
        false_alarms_method1 = false_alarms_method1 + 1;
    end

    % Método 2 (H0): Ruído Gaussiano com normalização std
    noise_H0_m2_raw = randn(total_samples, 1);
    noise_H0_m2 = noise_H0_m2_raw / std(noise_H0_m2_raw); % Normalização característica

    Y_H0_m2_all_windows = zeros(NFFT/2+1, M_windows);
    for k = 1:M_windows
        start_idx = (k-1)*NFFT + 1;
        end_idx = k*NFFT;
        if end_idx > total_samples; continue; end % Proteção
        window_data = noise_H0_m2(start_idx:end_idx);
        fft_window = fft(window_data.*hann(NFFT), NFFT);
        Y_H0_m2_all_windows(:, k) = fft_window(1:NFFT/2+1);
    end
    msc_H0_m2 = msc_fft(Y_H0_m2_all_windows, M_windows);
    if msc_H0_m2(sig_bin_idx) > threshold_msc
        false_alarms_method2 = false_alarms_method2 + 1;
    end
end
pfa_empirical_method1 = false_alarms_method1 / N_trials_H0;
pfa_empirical_method2 = false_alarms_method2 / N_trials_H0;

fprintf('PFA Empírica - Método 1: %.4f\n', pfa_empirical_method1);
fprintf('PFA Empírica - Método 2: %.4f\n', pfa_empirical_method2);

% --- Loop sobre diferentes SNRs para PD e MSC (H1: Sinal + Ruído) ---
pd_results_method1 = zeros(length(SNR_dB_range), 1);
pd_results_method2 = zeros(length(SNR_dB_range), 1);
msc_results_method1_H1 = zeros(length(SNR_dB_range), 1); % Para curva MSC vs SNR
msc_results_method2_H1 = zeros(length(SNR_dB_range), 1);

% Gerar Sinal Senoidal Puro uma vez
amplitude_signal_pure = 1;
signal_pure = amplitude_signal_pure * sin(2*pi*signal_freq*t);

fprintf('Calculando PD e MSC para H1 com %d trials por SNR...\n', N_trials_H1);
for snr_idx = 1:length(SNR_dB_range)
    current_SNR_dB = SNR_dB_range(snr_idx);
    current_SNR_linear = 10^(current_SNR_dB/10);

    detections_method1 = 0;
    detections_method2 = 0;
    temp_msc_sum_m1 = 0;
    temp_msc_sum_m2 = 0;

    for trial_idx = 1:N_trials_H1
        % --- Método 1: Simples (awgn) ---
        signal_noisy_method1 = awgn(signal_pure, current_SNR_dB, 'measured');
        % Opcional: Normalizar para std=1 para comparação mais direta de espectro
        signal_noisy_method1_norm = signal_noisy_method1 / std(signal_noisy_method1);
        signal_noisy_method1 = signal_noisy_method1_norm;
        Y_m1_all_windows = zeros(NFFT/2+1, M_windows);
        for k = 1:M_windows
            start_idx = (k-1)*NFFT + 1;
            end_idx = k*NFFT;
            window_data = signal_noisy_method1(start_idx:end_idx);
            fft_window = fft(window_data.*hann(NFFT), NFFT);
            Y_m1_all_windows(:, k) = fft_window(1:NFFT/2+1);
        end
        msc_values_m1 = msc_fft(Y_m1_all_windows, M_windows);
        current_msc_m1 = msc_values_m1(sig_bin_idx);
        temp_msc_sum_m1 = temp_msc_sum_m1 + current_msc_m1;
        if current_msc_m1 > threshold_msc
            detections_method1 = detections_method1 + 1;
        end
        if current_SNR_dB == SNR_dB_for_spectrum_plot && trial_idx == 1 && snr_idx == find(SNR_dB_range == SNR_dB_for_spectrum_plot,1)
            spectrum_method1_H1 = abs(Y_m1_all_windows(:,1)); % Espectro da primeira janela
        end

        % --- Método 2: Estilo rlord_gen_states ---
        power_signal_pure = var(signal_pure);
        if power_signal_pure == 0; power_signal_pure = eps; end
        power_noise_needed = power_signal_pure / current_SNR_linear;

        noise_raw_m2 = randn(size(signal_pure));
        noise_scaled_m2 = noise_raw_m2 * sqrt(power_noise_needed / var(noise_raw_m2));
        signal_noisy_method2 = signal_pure + noise_scaled_m2;
        signal_noisy_method2 = signal_noisy_method2 / std(signal_noisy_method2); % Normalização final

        Y_m2_all_windows = zeros(NFFT/2+1, M_windows);
        for k = 1:M_windows
            start_idx = (k-1)*NFFT + 1;
            end_idx = k*NFFT;
            window_data = signal_noisy_method2(start_idx:end_idx);
            fft_window = fft(window_data.*hann(NFFT), NFFT);
            Y_m2_all_windows(:, k) = fft_window(1:NFFT/2+1);
        end
        msc_values_m2 = msc_fft(Y_m2_all_windows, M_windows);
        current_msc_m2 = msc_values_m2(sig_bin_idx);
        temp_msc_sum_m2 = temp_msc_sum_m2 + current_msc_m2;
        if current_msc_m2 > threshold_msc
            detections_method2 = detections_method2 + 1;
        end
        if current_SNR_dB == SNR_dB_for_spectrum_plot && trial_idx == 1 && snr_idx == find(SNR_dB_range == SNR_dB_for_spectrum_plot,1)
             spectrum_method2_H1 = abs(Y_m2_all_windows(:,1)); % Espectro da primeira janela
        end
    end
    pd_results_method1(snr_idx) = detections_method1 / N_trials_H1;
    pd_results_method2(snr_idx) = detections_method2 / N_trials_H1;
    msc_results_method1_H1(snr_idx) = temp_msc_sum_m1 / N_trials_H1;
    msc_results_method2_H1(snr_idx) = temp_msc_sum_m2 / N_trials_H1;

    fprintf('Calculado para SNR = %.1f dB: PD_M1=%.3f, PD_M2=%.3f, MSC_M1=%.3f, MSC_M2=%.3f\n', ...
        current_SNR_dB, pd_results_method1(snr_idx), pd_results_method2(snr_idx), ...
        msc_results_method1_H1(snr_idx), msc_results_method2_H1(snr_idx));
end


% --- Plotar Resultados ---
% 1. Comparação dos Espectros de Frequência para H1 e uma SNR específica
if exist('spectrum_method1_H1', 'var') && exist('spectrum_method2_H1', 'var')
    figure;
    subplot(2,1,1);
    plot(f_axis, 20*log10(spectrum_method1_H1 + eps)); % Adicionado eps para evitar log10(0)
    title(['Método 1 (awgn): Espectro H1 para SNR = ' num2str(SNR_dB_for_spectrum_plot) ' dB']);
    xlabel('Frequência (Hz)'); ylabel('Magnitude (dB)'); grid on; xlim([0 FS/2]);
    yl = ylim;

    subplot(2,1,2);
    plot(f_axis, 20*log10(spectrum_method2_H1 + eps));
    title(['Método 2 (Estilo rlord): Espectro H1 para SNR = ' num2str(SNR_dB_for_spectrum_plot) ' dB']);
    xlabel('Frequência (Hz)'); ylabel('Magnitude (dB)'); grid on; xlim([0 FS/2]);
    ylim(yl);
    sgtitle(['Comparação dos Espectros de Frequência (Sinal + Ruído, F_{sig}=' num2str(signal_freq) ' Hz)']);
else
    disp('Espectros para plot não foram gerados. Verifique SNR_dB_for_spectrum_plot e SNR_dB_range.');
end

% 2. Curva MSC vs. SNR (para H1)
figure;
plot(SNR_dB_range, msc_results_method1_H1, 'b-o', 'DisplayName', 'Método 1 (awgn)');
hold on;
plot(SNR_dB_range, msc_results_method2_H1, 'r-s', 'DisplayName', 'Método 2 (Estilo rlord)');
hold off;
title(['MSC na Frequência do Sinal (' num2str(signal_freq) ' Hz) vs. SNR']);
xlabel('SNR (dB)'); ylabel('MSC Média'); legend show; grid on; ylim([0 1.05]);

% 3. Curva PD vs. SNR
figure;
plot(SNR_dB_range, pd_results_method1, 'b-o', 'DisplayName', 'PD - Método 1 (awgn)');
hold on;
plot(SNR_dB_range, pd_results_method2, 'r-s', 'DisplayName', 'PD - Método 2 (Estilo rlord)');
yline(1-target_PFA, 'k--', 'DisplayName', sprintf('1 - PFA Alvo = %.2f', 1-target_PFA));
yline(target_PFA, 'g--', 'DisplayName', sprintf('PFA Alvo = %.2f', target_PFA));
hold off;
title(['Probabilidade de Detecção (PD) vs. SNR (F_{sig}=' num2str(signal_freq) ' Hz, PFA Alvo=' num2str(target_PFA) ')']);
text(SNR_dB_range(1)+1, 0.15, sprintf('PFA Emp. M1: %.4f', pfa_empirical_method1), 'Color', 'blue');
text(SNR_dB_range(1)+1, 0.05, sprintf('PFA Emp. M2: %.4f', pfa_empirical_method2), 'Color', 'red');
xlabel('SNR (dB)'); ylabel('Probabilidade de Detecção (PD)');
legend show; grid on; ylim([-0.05 1.05]);

disp('Simulação e plots concluídos.');

% --- Função msc_fft (incluída para ser autocontida) ---
% function [ORD] = msc_fft(Y,M_win) % Renomeado M para M_win para clareza
% % Y é a FFT dos dados, com dimensão [NFFT/2+1, M_windows]
% % M_win é o número de janelas (M_windows)
% 
% if size(Y,2) == 0 % Sem janelas, sem MSC
%     ORD = zeros(size(Y,1),1);
%     return;
% end
% if M_win == 0
%     ORD = zeros(size(Y,1),1);
%     return;
% end
% 
% if (size(Y,2) ~= M_win) && M_win ~= 1
%     if M_win==1 && size(Y,2) == 1
%         ORD = ones(size(Y,1),1); % MSC é 1 para uma única janela
%         return;
%     else
%          % Se M_win for maior que o número de colunas em Y, isso é um erro.
%          % Mas se Y tiver menos colunas que M_win (ex: final do sinal),
%          % devemos usar o número real de janelas em Y.
%          if M_win > size(Y,2)
%              %warning('M_win (%d) é maior que o número de janelas em Y (%d). Usando size(Y,2).', M_win, size(Y,2));
%              M_win = size(Y,2);
%              if M_win == 0 % Se acabou sendo zero
%                  ORD = zeros(size(Y,1),1);
%                  return;
%              elseif M_win == 1
%                  ORD = ones(size(Y,1),1);
%                  return;
%              end
%          else
%             error('Número de janelas M_win inconsistente com as dimensões de Y.');
%          end
%     end
% end
% if M_win == 1 % Caso de uma única janela, MSC é 1.
%     ORD = ones(size(Y,1),1);
%     return;
% end
% 
% MSC
numerator = abs(sum(Y,2)).^2;
denominator = M_win * sum(abs(Y).^2,2);

ORD = numerator ./ denominator;
ORD(denominator == 0) = 0; % Se não há energia, MSC é 0 (ou NaN, mas 0 é mais prático aqui)
ORD(isnan(ORD)) = 0;
ORD(isinf(ORD)) = 0;
endx