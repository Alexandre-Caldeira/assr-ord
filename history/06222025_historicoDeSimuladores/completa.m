clear all;
close all;
clc;

% --- Parâmetros da Simulação ---
FS = 1000;             % Taxa de amostragem (Hz)
T_window = 0.5;        % Duração de cada janela para FFT/MSC (s)
NFFT = round(T_window * FS);  % Número de pontos da FFT (tamanho da janela)
M_windows = 10;        % Número de janelas para cálculo da MSC
signal_freq = 82;      % Frequência do sinal de interesse (Hz)

SNR_dB_range = -25:2:5; % Faixa de SNRs para testar (dB)
SNR_dB_for_spectrum_plot = -5; % dB, para plotar o espectro

% Parâmetros para PD e PFA
target_PFA = 0.05;     % PFA alvo para definir o limiar
N_trials_H0 = 10000;   % Número de realizações para estimar PFA empírica
N_trials_H1 = 2000;    % Número de realizações por SNR para estimar PD

% --- Calcular Limiar Teórico para MSC ---
threshold_msc = betainv(1 - target_PFA, 1, M_windows - 1);
fprintf('Limiar Teórico da MSC para PFA = %.3f: %.4f (M_windows=%d)\n', target_PFA, threshold_msc, M_windows);

% Duração total do sinal para M_windows para Métodos 1 e 2
total_duration = M_windows * T_window;
total_samples = round(total_duration * FS);
t = (0:total_samples-1)'/FS; % Vetor de tempo para o sinal completo

% Encontrar o bin de frequência mais próximo de signal_freq
f_axis_full = (0:NFFT/2)*FS/NFFT; % Eixo de frequência para FFT completa (0 a Nyquist)
[~, sig_bin_idx_full] = min(abs(f_axis_full - signal_freq));
fprintf('Frequência do sinal: %.1f Hz\n', signal_freq);
fprintf('  Bin FFT (completa, 0-idx para f_axis_full): %d (%.2f Hz)\n', sig_bin_idx_full-1, f_axis_full(sig_bin_idx_full));

% Ajuste do bin para o output de genSignals.m (que remove DC)
% genSignals retorna S_ com frequências Y(2:floor(end/2)+1) da FFT original.
% Se sig_bin_idx_full = 1 (DC), não estará em S_.
% Se sig_bin_idx_full > 1, o índice em S_ será sig_bin_idx_full - 1.
if sig_bin_idx_full == 1 && signal_freq ~= 0
    error('Frequência do sinal mapeada para DC, mas genSignals remove DC. Verifique signal_freq ou NFFT.');
end
sig_bin_idx_genSignals = sig_bin_idx_full - 1;
if sig_bin_idx_genSignals < 1 && signal_freq ~= 0 % Checagem adicional
    warning('Índice do sinal para genSignals é < 1. Verifique signal_freq (%.1f Hz) e NFFT (%d). sig_bin_idx_full=%d', signal_freq, NFFT, sig_bin_idx_full);
    % Pode ser necessário um tratamento mais robusto se signal_freq for muito baixo e NFFT pequeno
end
% O eixo de frequência para genSignals output (S3, S5)
% S_ tem floor(NFFT/2) pontos de frequência se NFFT é par.
% Corresponde a f_axis_full(2 : NFFT/2+1) se NFFT é par.
f_axis_genSignals = f_axis_full(2 : floor(NFFT/2)+1);
% Validação do sig_bin_idx_genSignals
if sig_bin_idx_genSignals > 0 && sig_bin_idx_genSignals <= length(f_axis_genSignals)
    fprintf('  Bin FFT (genSignals, 1-idx): %d (%.2f Hz)\n', sig_bin_idx_genSignals, f_axis_genSignals(sig_bin_idx_genSignals));
else
    if signal_freq == 0 % Se o sinal for DC, ele não estará no output de genSignals
        fprintf('  Aviso: Frequência do sinal é DC (0 Hz), que é removido por genSignals.\n');
        % Não poderemos calcular MSC para o Método 3 se o sinal for DC.
    else
        error('Índice do sinal sig_bin_idx_genSignals (%d) fora do intervalo para f_axis_genSignals (1 a %d).', sig_bin_idx_genSignals, length(f_axis_genSignals));
    end
end


% --- Estimar PFA Empírica (H0: Somente Ruído) ---
false_alarms_method1 = 0;
false_alarms_method2 = 0;
false_alarms_method3 = 0;

fprintf('Estimando PFA empírica com %d trials...\n', N_trials_H0);

% Método 3 (genSignals) para H0
% SNRfun para H0 - o valor exato da SNR não deve importar para S3, pois é ruído puro
snr_dummy_for_H0_genSignals = @() -100; % SNR em dB muito baixa
% genSignals retorna S_(:,:,ii) como [Nfreqs, Njanelas]
% Nsinais em genSignals é N_trials_H0
[~, ~, S3_all_trials, ~, ~, ~] = genSignals(snr_dummy_for_H0_genSignals, FS, signal_freq, NFFT, N_trials_H0, M_windows);

for trial_idx = 1:N_trials_H0
    % Método 1 (H0)
    noise_H0_m1 = randn(total_samples, 1);
    Y_H0_m1_all_windows = zeros(NFFT/2+1, M_windows);
    for k = 1:M_windows
        start_idx = (k-1)*NFFT + 1; end_idx = k*NFFT;
        if end_idx > total_samples; continue; end
        window_data = noise_H0_m1(start_idx:end_idx);
        fft_window = fft(window_data.*hann(NFFT), NFFT);
        Y_H0_m1_all_windows(:, k) = fft_window(1:NFFT/2+1);
    end
    msc_H0_m1 = msc_fft(Y_H0_m1_all_windows, M_windows);
    if ~isempty(msc_H0_m1) && msc_H0_m1(sig_bin_idx_full) > threshold_msc
        false_alarms_method1 = false_alarms_method1 + 1;
    end

    % Método 2 (H0)
    noise_H0_m2_raw = randn(total_samples, 1);
    noise_H0_m2 = noise_H0_m2_raw / std(noise_H0_m2_raw);
    Y_H0_m2_all_windows = zeros(NFFT/2+1, M_windows);
    for k = 1:M_windows
        start_idx = (k-1)*NFFT + 1; end_idx = k*NFFT;
        if end_idx > total_samples; continue; end
        window_data = noise_H0_m2(start_idx:end_idx);
        fft_window = fft(window_data.*hann(NFFT), NFFT);
        Y_H0_m2_all_windows(:, k) = fft_window(1:NFFT/2+1);
    end
    msc_H0_m2 = msc_fft(Y_H0_m2_all_windows, M_windows);
     if ~isempty(msc_H0_m2) && msc_H0_m2(sig_bin_idx_full) > threshold_msc
        false_alarms_method2 = false_alarms_method2 + 1;
    end

    % Método 3 (H0) - usando S3 de genSignals
    if signal_freq ~= 0 % Só calcular se o sinal não for DC
        Y_H0_m3_single_trial = S3_all_trials(:,:,trial_idx); % [Nfreqs_genSignals, M_windows]
        msc_H0_m3 = msc_fft(Y_H0_m3_single_trial, M_windows);
        if ~isempty(msc_H0_m3) && msc_H0_m3(sig_bin_idx_genSignals) > threshold_msc
            false_alarms_method3 = false_alarms_method3 + 1;
        end
    end
end
pfa_empirical_method1 = false_alarms_method1 / N_trials_H0;
pfa_empirical_method2 = false_alarms_method2 / N_trials_H0;
if signal_freq ~= 0
    pfa_empirical_method3 = false_alarms_method3 / N_trials_H0;
    fprintf('PFA Empírica - Método 3 (genSignals): %.4f\n', pfa_empirical_method3);
else
    pfa_empirical_method3 = NaN;
    fprintf('PFA Empírica - Método 3 (genSignals): N/A (sinal em DC)\n');
end
fprintf('PFA Empírica - Método 1: %.4f\n', pfa_empirical_method1);
fprintf('PFA Empírica - Método 2: %.4f\n', pfa_empirical_method2);


% --- Loop sobre diferentes SNRs para PD e MSC (H1: Sinal + Ruído) ---
pd_results_method1 = zeros(length(SNR_dB_range), 1);
pd_results_method2 = zeros(length(SNR_dB_range), 1);
pd_results_method3 = zeros(length(SNR_dB_range), 1);
msc_results_method1_H1 = zeros(length(SNR_dB_range), 1);
msc_results_method2_H1 = zeros(length(SNR_dB_range), 1);
msc_results_method3_H1 = zeros(length(SNR_dB_range), 1);

% Gerar Sinal Senoidal Puro uma vez para Métodos 1 e 2
amplitude_signal_pure = 1;
signal_pure_m1_m2 = amplitude_signal_pure * sin(2*pi*signal_freq*t);

fprintf('Calculando PD e MSC para H1 com %d trials por SNR...\n', N_trials_H1);
for snr_idx = 1:length(SNR_dB_range)
    current_SNR_dB = SNR_dB_range(snr_idx);
    current_SNR_linear = 10^(current_SNR_dB/10);

    detections_method1 = 0; temp_msc_sum_m1 = 0;
    detections_method2 = 0; temp_msc_sum_m2 = 0;
    detections_method3 = 0; temp_msc_sum_m3 = 0;

    % Método 3 (genSignals) para H1
    if signal_freq ~=0
        SNRfun_H1_genSignals = @() current_SNR_dB; % genSignals espera SNR em dB
        [~,~,~,~,S5_all_trials_m3,~] = genSignals(SNRfun_H1_genSignals, FS, signal_freq, NFFT, N_trials_H1, M_windows);
    end

    for trial_idx = 1:N_trials_H1
        % Método 1 (awgn)
        signal_noisy_m1 = awgn(signal_pure_m1_m2, current_SNR_dB, 'measured');
        Y_m1_all_windows = zeros(NFFT/2+1, M_windows);
        for k = 1:M_windows
            start_idx = (k-1)*NFFT + 1; end_idx = k*NFFT;
            window_data = signal_noisy_m1(start_idx:end_idx);
            fft_window = fft(window_data.*hann(NFFT), NFFT);
            Y_m1_all_windows(:, k) = fft_window(1:NFFT/2+1);
        end
        msc_values_m1 = msc_fft(Y_m1_all_windows, M_windows);
        if ~isempty(msc_values_m1)
            current_msc_m1 = msc_values_m1(sig_bin_idx_full);
            temp_msc_sum_m1 = temp_msc_sum_m1 + current_msc_m1;
            if current_msc_m1 > threshold_msc; detections_method1 = detections_method1 + 1; end
        end
        if current_SNR_dB == SNR_dB_for_spectrum_plot && trial_idx == 1 && snr_idx == find(SNR_dB_range == SNR_dB_for_spectrum_plot,1) && ~isempty(msc_values_m1)
            spectrum_method1_H1 = abs(Y_m1_all_windows(:,1));
        end

        % Método 2 (Estilo rlord_gen_states)
        power_signal_pure = var(signal_pure_m1_m2);
        if power_signal_pure == 0; power_signal_pure = eps; end
        power_noise_needed = power_signal_pure / current_SNR_linear;
        noise_raw_m2 = randn(size(signal_pure_m1_m2));
        noise_scaled_m2 = noise_raw_m2 * sqrt(power_noise_needed / var(noise_raw_m2));
        signal_noisy_m2 = signal_pure_m1_m2 + noise_scaled_m2;
        signal_noisy_m2 = signal_noisy_m2 / std(signal_noisy_m2);
        Y_m2_all_windows = zeros(NFFT/2+1, M_windows);
        for k = 1:M_windows
            start_idx = (k-1)*NFFT + 1; end_idx = k*NFFT;
            window_data = signal_noisy_m2(start_idx:end_idx);
            fft_window = fft(window_data.*hann(NFFT), NFFT);
            Y_m2_all_windows(:, k) = fft_window(1:NFFT/2+1);
        end
        msc_values_m2 = msc_fft(Y_m2_all_windows, M_windows);
        if ~isempty(msc_values_m2)
            current_msc_m2 = msc_values_m2(sig_bin_idx_full);
            temp_msc_sum_m2 = temp_msc_sum_m2 + current_msc_m2;
            if current_msc_m2 > threshold_msc; detections_method2 = detections_method2 + 1; end
        end
        if current_SNR_dB == SNR_dB_for_spectrum_plot && trial_idx == 1 && snr_idx == find(SNR_dB_range == SNR_dB_for_spectrum_plot,1) && ~isempty(msc_values_m2)
             spectrum_method2_H1 = abs(Y_m2_all_windows(:,1));
        end

        % Método 3 (genSignals)
        if signal_freq ~=0
            Y_m3_single_trial = S5_all_trials_m3(:,:,trial_idx);
            msc_values_m3 = msc_fft(Y_m3_single_trial, M_windows);
            if ~isempty(msc_values_m3)
                current_msc_m3 = msc_values_m3(sig_bin_idx_genSignals);
                temp_msc_sum_m3 = temp_msc_sum_m3 + current_msc_m3;
                if current_msc_m3 > threshold_msc; detections_method3 = detections_method3 + 1; end
            end
            if current_SNR_dB == SNR_dB_for_spectrum_plot && trial_idx == 1 && snr_idx == find(SNR_dB_range == SNR_dB_for_spectrum_plot,1) && ~isempty(msc_values_m3)
                 spectrum_method3_H1_genSignalBins = abs(Y_m3_single_trial(:,1));
            end
        end
    end
    pd_results_method1(snr_idx) = detections_method1 / N_trials_H1;
    pd_results_method2(snr_idx) = detections_method2 / N_trials_H1;
    msc_results_method1_H1(snr_idx) = temp_msc_sum_m1 / N_trials_H1;
    msc_results_method2_H1(snr_idx) = temp_msc_sum_m2 / N_trials_H1;
    if signal_freq ~= 0
        pd_results_method3(snr_idx) = detections_method3 / N_trials_H1;
        msc_results_method3_H1(snr_idx) = temp_msc_sum_m3 / N_trials_H1;
        fprintf('SNR=%.1fdB: PD(M1)=%.3f,PD(M2)=%.3f,PD(M3)=%.3f | MSC(M1)=%.3f,MSC(M2)=%.3f,MSC(M3)=%.3f\n', ...
            current_SNR_dB, pd_results_method1(snr_idx), pd_results_method2(snr_idx), pd_results_method3(snr_idx), ...
            msc_results_method1_H1(snr_idx), msc_results_method2_H1(snr_idx), msc_results_method3_H1(snr_idx));
    else
        pd_results_method3(snr_idx) = NaN;
        msc_results_method3_H1(snr_idx) = NaN;
        fprintf('SNR=%.1fdB: PD(M1)=%.3f,PD(M2)=%.3f,PD(M3)=N/A | MSC(M1)=%.3f,MSC(M2)=%.3f,MSC(M3)=N/A\n', ...
            current_SNR_dB, pd_results_method1(snr_idx), pd_results_method2(snr_idx), ...
            msc_results_method1_H1(snr_idx), msc_results_method2_H1(snr_idx));
    end
end


% --- Plotar Resultados ---
% 1. Espectros de Frequência para H1
figure_spectra = figure;
num_methods_for_spectra = 2;
if exist('spectrum_method3_H1_genSignalBins', 'var'); num_methods_for_spectra = 3; end

subplot_idx = 1;
if exist('spectrum_method1_H1', 'var')
    subplot(num_methods_for_spectra,1,subplot_idx); subplot_idx = subplot_idx + 1;
    plot(f_axis_full, 20*log10(spectrum_method1_H1 + eps));
    title(['M1 (awgn): Espectro H1, SNR = ' num2str(SNR_dB_for_spectrum_plot) ' dB']);
    xlabel('Frequência (Hz)'); ylabel('Magnitude (dB)'); grid on; xlim([0 FS/2]);
    yl_spec = ylim;
end
if exist('spectrum_method2_H1', 'var')
    subplot(num_methods_for_spectra,1,subplot_idx); subplot_idx = subplot_idx + 1;
    plot(f_axis_full, 20*log10(spectrum_method2_H1 + eps));
    title(['M2 (Estilo rlord): Espectro H1, SNR = ' num2str(SNR_dB_for_spectrum_plot) ' dB']);
    xlabel('Frequência (Hz)'); ylabel('Magnitude (dB)'); grid on; xlim([0 FS/2]);
    if exist('yl_spec', 'var'); ylim(yl_spec); end
end
if exist('spectrum_method3_H1_genSignalBins', 'var') && signal_freq ~= 0
    subplot(num_methods_for_spectra,1,subplot_idx);
    plot(f_axis_genSignals, 20*log10(spectrum_method3_H1_genSignalBins + eps)); % Usa f_axis_genSignals
    title(['M3 (genSignals): Espectro H1, SNR = ' num2str(SNR_dB_for_spectrum_plot) ' dB']);
    xlabel('Frequência (Hz)'); ylabel('Magnitude (dB)'); grid on; xlim([0 FS/2]);
     if exist('yl_spec', 'var'); ylim(yl_spec); end
end
if subplot_idx > 1
    sgtitle(['Comparação Espectros (Sinal + Ruído, F_{sig}=' num2str(signal_freq) ' Hz)']);
else
    disp('Espectros para plot não foram gerados.');
    close(figure_spectra); % Fecha a figura se não houver nada para plotar
end


% 2. Curva MSC vs. SNR (para H1)
figure;
plot(SNR_dB_range, msc_results_method1_H1, 'b-o', 'DisplayName', 'M1 (awgn)');
hold on;
plot(SNR_dB_range, msc_results_method2_H1, 'r-s', 'DisplayName', 'M2 (Estilo rlord)');
if signal_freq ~= 0
    plot(SNR_dB_range, msc_results_method3_H1, 'g-^', 'DisplayName', 'M3 (genSignals)');
end
hold off;
title(['MSC na Frequência do Sinal (' num2str(signal_freq) ' Hz) vs. SNR']);
xlabel('SNR (dB)'); ylabel('MSC Média'); legend show; grid on; ylim([0 1.05]);

% 3. Curva PD vs. SNR
figure;
plot(SNR_dB_range, pd_results_method1, 'b-o', 'DisplayName', 'PD - M1 (awgn)');
hold on;
plot(SNR_dB_range, pd_results_method2, 'r-s', 'DisplayName', 'PD - M2 (Estilo rlord)');
if signal_freq ~= 0
    plot(SNR_dB_range, pd_results_method3, 'g-^', 'DisplayName', 'PD - M3 (genSignals)');
end
yline(1-target_PFA, 'k--', 'DisplayName', sprintf('1 - PFA Alvo = %.2f', 1-target_PFA));
yline(target_PFA, 'g--', 'DisplayName', sprintf('PFA Alvo = %.2f', target_PFA)); % Esta linha pode confundir com a curva PD
hold off;
title_str_pd = sprintf('PD vs. SNR (F_{sig}=%.1f Hz, PFA Alvo=%.2f)', signal_freq, target_PFA);
texts_pfa = {sprintf('PFA Emp. M1: %.4f', pfa_empirical_method1), ...
             sprintf('PFA Emp. M2: %.4f', pfa_empirical_method2)};
if signal_freq ~=0
    texts_pfa{end+1} = sprintf('PFA Emp. M3: %.4f', pfa_empirical_method3);
end
% Posicionar textos de PFA
text_y_start = 0.25; text_y_step = -0.07;
for i=1:length(texts_pfa)
    text(SNR_dB_range(1)+0.5, text_y_start + (i-1)*text_y_step, texts_pfa{i});
end
title(title_str_pd);
xlabel('SNR (dB)'); ylabel('Probabilidade de Detecção (PD)');
legend show; grid on; ylim([-0.05 1.05]);

disp('Simulação e plots concluídos.');

% --- genSignals.m (COLOQUE EM UM ARQUIVO SEPARADO 'genSignals.m') ---
% function [S1, S2, S3, S4, S5, S6] = genSignals(SNRfun, FS, SFREQ, NFFT, Nsinais, Njanelas)
% % ... (conteúdo do seu genSignals.m aqui) ...
% end

% --- Função msc_fft (incluída para ser autocontida, ou coloque em arquivo separado) ---
% (Função msc_fft como definida no início desta resposta)
% function [ORD] = msc_fft(Y_input, M_win) ... end