clear all;
close all;
clc;

% --- Parâmetros da Simulação ---
FS = 1000;             % Taxa de amostragem (Hz)
T_window = 0.5;        % Duração de cada janela para FFT/MSC (s)
NFFT = T_window * FS;  % Número de pontos da FFT (igual ao tamanho da janela)
M_windows = 20;        % Número de janelas para cálculo da MSC
signal_freq = 50;      % Frequência do sinal de interesse (Hz)

SNR_dB_range = -25:2:10; % Faixa de SNRs para testar (dB)
num_realizations = 50;   % Número de realizações para promediar a MSC (para suavizar a curva)

% Para plotar o espectro, escolha uma SNR
SNR_dB_for_spectrum_plot = -5; % dB

% Inicializar arrays para armazenar resultados da MSC
msc_results_method1 = zeros(length(SNR_dB_range), 1);
msc_results_method2 = zeros(length(SNR_dB_range), 1);

% --- Loop sobre diferentes SNRs ---
for snr_idx = 1:length(SNR_dB_range)
    current_SNR_dB = SNR_dB_range(snr_idx);
    current_SNR_linear = 10^(current_SNR_dB/10);

    temp_msc_method1 = zeros(num_realizations, 1);
    temp_msc_method2 = zeros(num_realizations, 1);

    for real_idx = 1:num_realizations
        % Duração total do sinal para M_windows
        total_duration = M_windows * T_window;
        t = (0:1/FS:total_duration-1/FS)'; % Vetor de tempo para o sinal completo

        % 1. Gerar Sinal Senoidal Puro
        amplitude_signal_pure = 1; % Amplitude arbitrária para o sinal puro inicial
        signal_pure = amplitude_signal_pure * sin(2*pi*signal_freq*t);

        % --- Método 1: Simples (awgn) ---
        signal_noisy_method1 = awgn(signal_pure, current_SNR_dB, 'measured');
        % Opcional: normalizar como em rlord_gen_states (para comparação mais direta de espectros)
        % signal_noisy_method1 = signal_noisy_method1 / std(signal_noisy_method1);

        % Segmentar e calcular FFT para Método 1
        Y_method1_all_windows = zeros(NFFT/2+1, M_windows);
        for k = 1:M_windows
            start_idx = (k-1)*NFFT + 1;
            end_idx = k*NFFT;
            window_data = signal_noisy_method1(start_idx:end_idx);
            fft_window = fft(window_data.*hann(NFFT), NFFT); % Janela de Hann para reduzir vazamento
            Y_method1_all_windows(:, k) = fft_window(1:NFFT/2+1);
        end

        % Calcular MSC para Método 1
        msc_values_method1 = msc_fft(Y_method1_all_windows, M_windows);
        f_axis = (0:NFFT/2)*FS/NFFT;
        [~, sig_bin_idx] = min(abs(f_axis - signal_freq));
        temp_msc_method1(real_idx) = msc_values_method1(sig_bin_idx);

        % Guardar espectro para uma realização específica (para o plot de espectro)
        if current_SNR_dB == SNR_dB_for_spectrum_plot && real_idx == 1
            spectrum_method1 = abs(Y_method1_all_windows(:,1)); % Espectro da primeira janela
        end


        % --- Método 2: Estilo rlord_gen_states ---
        % Potência do sinal puro
        power_signal_pure = var(signal_pure); % ou (amplitude_signal_pure^2)/2 se não houver DC
        if power_signal_pure == 0; power_signal_pure = eps; end % Evitar divisão por zero se sinal for nulo

        % Potência de ruído necessária
        power_noise_needed = power_signal_pure / current_SNR_linear;

        % Gerar ruído gaussiano e escalar
        noise_raw = randn(size(signal_pure));
        noise_scaled = noise_raw * sqrt(power_noise_needed / var(noise_raw));

        signal_noisy_method2 = signal_pure + noise_scaled;
        % Normalização final como em rlord_gen_states
        signal_noisy_method2 = signal_noisy_method2 / std(signal_noisy_method2);


        % Segmentar e calcular FFT para Método 2
        Y_method2_all_windows = zeros(NFFT/2+1, M_windows);
        for k = 1:M_windows
            start_idx = (k-1)*NFFT + 1;
            end_idx = k*NFFT;
            window_data = signal_noisy_method2(start_idx:end_idx);
            fft_window = fft(window_data.*hann(NFFT), NFFT); % Janela de Hann
            Y_method2_all_windows(:, k) = fft_window(1:NFFT/2+1);
        end

        % Calcular MSC para Método 2
        msc_values_method2 = msc_fft(Y_method2_all_windows, M_windows);
        temp_msc_method2(real_idx) = msc_values_method2(sig_bin_idx);

        % Guardar espectro para uma realização específica (para o plot de espectro)
        if current_SNR_dB == SNR_dB_for_spectrum_plot && real_idx == 1
            spectrum_method2 = abs(Y_method2_all_windows(:,1)); % Espectro da primeira janela
        end
    end

    msc_results_method1(snr_idx) = mean(temp_msc_method1);
    msc_results_method2(snr_idx) = mean(temp_msc_method2);

    fprintf('Calculado para SNR = %.1f dB\n', current_SNR_dB);
end

% --- Plotar Resultados ---

% 1. Comparação dos Espectros de Frequência para uma SNR específica
figure;
subplot(2,1,1);
plot(f_axis, 20*log10(spectrum_method1));
title(['Método 1 (awgn): Espectro para SNR = ' num2str(SNR_dB_for_spectrum_plot) ' dB']);
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
grid on;
xlim([0 FS/2]);
yl = ylim; % Para usar o mesmo limite y no próximo subplot

subplot(2,1,2);
plot(f_axis, 20*log10(spectrum_method2));
title(['Método 2 (Estilo rlord\_gen\_states): Espectro para SNR = ' num2str(SNR_dB_for_spectrum_plot) ' dB']);
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
grid on;
xlim([0 FS/2]);
ylim(yl); % Garante mesma escala Y para comparação visual

sgtitle('Comparação dos Espectros de Frequência (Primeira Janela)');

% 2. Curva MSC vs. SNR
figure;
plot(SNR_dB_range, msc_results_method1, 'b-o', 'DisplayName', 'Método 1 (awgn)');
hold on;
plot(SNR_dB_range, msc_results_method2, 'r-s', 'DisplayName', 'Método 2 (Estilo rlord\_gen\_states)');
hold off;
title(['MSC na Frequência do Sinal (' num2str(signal_freq) ' Hz) vs. SNR']);
xlabel('SNR (dB)');
ylabel('MSC');
legend show;
grid on;
ylim([0 1.05]);

disp('Simulação e plots concluídos.');

% --- START OF msc_fft.m (copiado aqui para ser autocontido) ---
function [ORD] = msc_fft(Y,M)
if (size(Y,2) ~= M) && M ~= 1
    if M==1
        ORD = ones(size(Y,1),1);
        return;
    else
        error('Número de janelas M inconsistente com as dimensões de Y.');
    end
end
ORD =  abs(sum(Y,2)).^2 ./ (M*sum(abs(Y).^2,2));
ORD(isnan(ORD)) = 0;
ORD(isinf(ORD)) = 0; % Adicionado para cobrir casos onde o denominador é exatamente zero
end
% --- END OF msc_fft.m ---

% --- START OF msc_fft.m (ou coloque no mesmo arquivo do script principal se preferir) ---
% function [ORD] = msc_fft(Y,M)
% % Y é a FFT dos dados, com dimensão [NFFT/2+1, M_windows]
% % M é o número de janelas (M_windows)
% 
% if (size(Y,2) ~= M) && M ~= 1 % M=1 é um caso especial, mas MSC não é bem definido
%     if M==1 % se for apenas uma janela, MSC é trivialmente 1 (ou indefinido)
%         % warning('MSC com M=1. Retornando 1 para todas as frequências.');
%         ORD = ones(size(Y,1),1);
%         return;
%     else
%         error('Número de janelas M inconsistente com as dimensões de Y.');
%     end
% end

% MSC
% abs(sum(Y,2)).^2 é o quadrado da magnitude da soma coerente das FFTs através das janelas
% M*sum(abs(Y).^2,2) é M vezes a soma das potências (magnitude ao quadrado) através das janelas
% ORD =  abs(sum(Y,2)).^2 ./ (M*sum(abs(Y).^2,2));
% % Lidar com possível divisão por zero se o denominador for zero (sinal zero em todas as janelas)
% ORD(isnan(ORD)) = 0; % Ou algum outro valor apropriado, como NaN
% end
% --- END OF msc_fft.m ---