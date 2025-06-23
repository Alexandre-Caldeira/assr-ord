clear all; close all; clc;

% --- Parâmetros ---
FS = 1000;                  % Taxa de amostragem (Hz)
NFFT = 1024;                % Número de pontos da FFT (e duração da janela de análise)
duration = NFFT / FS;       % Duração do sinal em segundos
t = (0:NFFT-1)' / FS;       % Vetor de tempo

f_signal = 100;             % Frequência do sinal (Hz)
% Para garantir que o sinal caia em um bin exato (opcional, mas simplifica a análise)
signal_bin_ideal_idx = round(f_signal * NFFT / FS);
f_signal = signal_bin_ideal_idx * FS / NFFT;
fprintf('Frequência do sinal ajustada para: %.2f Hz (bin %d)\n', f_signal, signal_bin_ideal_idx);


target_snr_db_per_bin = 10; % SNR desejada no bin do sinal (dB)
target_snr_linear_per_bin = 10^(target_snr_db_per_bin / 10);

% --- Geração de Sinal e Ruído ---

% 1. Fixar variância do ruído no tempo
sigma_n_sq = 1; % Variância do ruído no tempo (potência do ruído no tempo)
std_n = sqrt(sigma_n_sq);

% 2. Potência de ruído esperada por bin (não DC, não Nyquist) na FFT
% P_n_bin_expected = 2 * sigma_n_sq / NFFT; % Para um espectro unilateral
% Derivação alternativa: |fft(ruido_unit_var)(k)|^2 ~ NFFT. Potência no bin = (2/NFFT^2) * |fft(ruido_unit_var)(k)|^2 ~ 2/NFFT
P_n_bin_expected = (2/NFFT) * sigma_n_sq; % Potência em um bin da FFT para ruído com sigma_n_sq no tempo

% 3. Potência de sinal necessária no bin do sinal
P_s_bin_target = target_snr_linear_per_bin * P_n_bin_expected;

% 4. Amplitude A da senoide para atingir P_s_bin_target
% P_s_bin_target = A^2 / 2  (assumindo que toda a potência da senoide cai em um bin)
A_signal = sqrt(2 * P_s_bin_target);

% 5. Gerar componentes
s_t = A_signal * sin(2 * pi * f_signal * t);
n_t = std_n * randn(NFFT, 1); % Ruído gaussiano com std_n
x_t = s_t + n_t;              % Sinal com ruído

% --- Análise e Verificação ---
X_k = fft(x_t);
S_k = fft(s_t); % FFT do sinal puro para referência
N_k = fft(n_t); % FFT do ruído puro para referência

% Índices dos bins da FFT (0 a NFFT-1)
% Frequências para o espectro unilateral (excluindo DC para P_n_bin_expected)
freq_axis = (0:NFFT/2) * FS / NFFT; % Para plotagem
idx_dc = 1;
idx_nyquist = NFFT/2 + 1;

% Potência nos bins (usando a fórmula (2/NFFT^2)*|X(k)|^2 para unilateral)
power_X_k_onesided = (2/NFFT^2) * abs(X_k(1:idx_nyquist)).^2;
power_X_k_onesided(idx_dc) = (1/NFFT^2) * abs(X_k(idx_dc)).^2; % DC não é multiplicado por 2
if mod(NFFT,2) == 0 % Se NFFT é par, Nyquist também não é multiplicado por 2
    power_X_k_onesided(idx_nyquist) = (1/NFFT^2) * abs(X_k(idx_nyquist)).^2;
end

% Encontrar o bin correspondente à f_signal (deveria ser signal_bin_ideal_idx + 1 para array MATLAB)
actual_signal_bin_fft_idx = signal_bin_ideal_idx + 1; % MATLAB é 1-indexado, FFT é 0-indexado na teoria

measured_signal_power_in_bin = power_X_k_onesided(actual_signal_bin_fft_idx);

% Estimar potência de ruído nos bins adjacentes (ou média de vários bins de ruído)
noise_bins_for_avg = [actual_signal_bin_fft_idx+5, actual_signal_bin_fft_idx+6, actual_signal_bin_fft_idx-5, actual_signal_bin_fft_idx-6];
noise_bins_for_avg = noise_bins_for_avg(noise_bins_for_avg > idx_dc & noise_bins_for_avg < idx_nyquist & noise_bins_for_avg ~= actual_signal_bin_fft_idx); % Validar índices
if isempty(noise_bins_for_avg)
    noise_bins_for_avg = actual_signal_bin_fft_idx+2; % Fallback simples
    if noise_bins_for_avg >= idx_nyquist; noise_bins_for_avg = actual_signal_bin_fft_idx-2; end
end
measured_noise_power_around_signal = mean(power_X_k_onesided(noise_bins_for_avg));

% SNR por bin medida
measured_snr_per_bin = measured_signal_power_in_bin / measured_noise_power_around_signal;
measured_snr_db_per_bin = 10 * log10(measured_snr_per_bin);

fprintf('--- Verificação ---\n');
fprintf('Amplitude A do sinal: %.4f\n', A_signal);
fprintf('Potência esperada do sinal no bin (A^2/2): %.2e\n', A_signal^2/2);
fprintf('Potência esperada do ruído por bin (2*sigma_n^2/NFFT): %.2e\n', P_n_bin_expected);
fprintf('SNR por bin alvo: %.2f dB (Linear: %.2f)\n', target_snr_db_per_bin, target_snr_linear_per_bin);
fprintf('Potência medida do sinal no bin %d (freq %.2f Hz): %.2e\n', actual_signal_bin_fft_idx, freq_axis(actual_signal_bin_fft_idx), measured_signal_power_in_bin);
fprintf('Potência média medida do ruído em bins adjacentes: %.2e\n', measured_noise_power_around_signal);
fprintf('SNR por bin medida: %.2f dB (Linear: %.2f)\n', measured_snr_db_per_bin, measured_snr_per_bin);

% Plot
figure;
plot(freq_axis, 10*log10(power_X_k_onesided + eps)); % Adiciona eps para evitar log10(0)
hold on;
plot(freq_axis(actual_signal_bin_fft_idx), 10*log10(measured_signal_power_in_bin + eps), 'ro', 'MarkerSize', 8, 'DisplayName', 'Sinal no Bin');
plot(freq_axis(noise_bins_for_avg), 10*log10(power_X_k_onesided(noise_bins_for_avg) + eps), 'gx', 'MarkerSize', 8, 'DisplayName', 'Ruído Bins');
title(sprintf('Espectro de Potência Unilateral (SNR alvo por bin: %.1f dB)', target_snr_db_per_bin));
xlabel('Frequência (Hz)');
ylabel('Potência/Hz (dB)');
legend;
grid on;
xlim([0 FS/2]);