import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, scoreatpercentile, beta, chi2
import time

# --- Helper Functions ---
def msc_fft_py(Y_input, M_win):
    """
    Calcula a Magnitude Squared Coherence (MSC).
    Y_input: array 2D de FFTs (Nfreqs x M_windows)
    M_win: número de janelas
    """
    if Y_input.shape[1] == 0 or M_win == 0:
        return np.zeros(Y_input.shape[0])

    current_M_win = min(M_win, Y_input.shape[1])
    if current_M_win == 0:
        return np.zeros(Y_input.shape[0])

    if current_M_win == 1:
        return np.ones(Y_input.shape[0])

    Y_to_use = Y_input[:, :current_M_win]
    
    numerator = np.abs(np.sum(Y_to_use, axis=1))**2
    denominator = current_M_win * np.sum(np.abs(Y_to_use)**2, axis=1)
    
    msc = np.zeros_like(numerator, dtype=float)
    valid_den = denominator != 0
    msc[valid_den] = numerator[valid_den] / denominator[valid_den]
    msc[np.isnan(msc)] = 0
    msc[np.isinf(msc)] = 0
    return msc

def mcsm_py(y, tj, fs):
    """
    Calcula a Multiple Component Synchrony Measure (MCSM).
    y: sinal no domínio do tempo [total_samples, N_channels]
    tj: número de pontos por época (NFFT)
    fs: taxa de amostragem (não usado no cálculo principal)
    """
    tamsinal, N_channels = y.shape
    nfft = int(tj / 2)
    M_windows = int(tamsinal / tj)
    
    y = y[:M_windows * tj, :] # Garante que y seja um múltiplo exato de tj*M_windows

    # Remodela y para [tj, M_windows, N_channels] e aplica FFT ao longo do primeiro eixo (tj)
    Y_fft_per_window_channel = np.fft.fft(y.reshape(tj, M_windows, N_channels), axis=0)
    Y_fft_per_window_channel = Y_fft_per_window_channel[:nfft + 1, :, :] # [Nfreqs, M_windows, N_channels]

    teta = np.angle(Y_fft_per_window_channel)
    C = np.cos(teta)
    S = np.sin(teta)

    # Calcula a média de C e S ao longo dos canais (eixo=2)
    Cmed = np.mean(C, axis=2) # [Nfreqs, M_windows]
    Smed = np.mean(S, axis=2) # [Nfreqs, M_windows]

    teta_med = np.arctan2(Smed, Cmed) # [Nfreqs, M_windows]

    # Calcula MCSM_k = (1/M^2) * ( (sum_m cos(theta_k,m))^2 + (sum_m sin(theta_k,m))^2 )
    # Soma sobre as janelas (eixo=1)
    csmN = (1 / M_windows**2) * (np.sum(np.cos(teta_med), axis=1)**2 + np.sum(np.sin(teta_med), axis=1)**2)
    
    # Define os componentes DC e Nyquist como NaN, conforme o código MATLAB
    csmN[0] = np.nan # Componente DC
    if tj % 2 == 0: # Se NFFT for par, Nyquist está em nfft
        csmN[nfft] = np.nan # Componente Nyquist

    return csmN

def mgbt_py(y, x, tj, fs):
    """
    Calcula o Multivariate Global Beta Test (MGBT).
    y: sinal [total_samples, N_channels]
    x: ruído [total_samples, N_channels]
    tj: número de pontos por época (NFFT)
    fs: taxa de amostragem (não usado no cálculo principal)
    """
    nf = int(tj / 2) + 1
    N_channels_y = y.shape[1]
    N_channels_x = x.shape[1]

    if N_channels_y != N_channels_x:
        raise ValueError('O número de canais no sinal (y) e no ruído (x) deve ser o mesmo.')
    N_channels = N_channels_y

    M_windows_y = int(y.shape[0] / tj)
    M_windows_x = int(x.shape[0] / tj)

    y = y[:M_windows_y * tj, :]
    x = x[:M_windows_x * tj, :]

    y_reshaped = y.reshape(tj, M_windows_y, N_channels)
    x_reshaped = x.reshape(tj, M_windows_x, N_channels)

    # Normalização pelo desvio padrão para cada janela e canal
    std_y_per_window = np.std(y_reshaped, axis=0, keepdims=True)
    std_x_per_window = np.std(x_reshaped, axis=0, keepdims=True)

    std_y_per_window[std_y_per_window == 0] = 1
    std_x_per_window[std_x_per_window == 0] = 1

    y_normalized = y_reshaped / std_y_per_window
    x_normalized = x_reshaped / std_x_per_window

    Y_fft_sq = np.abs(np.fft.fft(y_normalized, axis=0))**2
    X_fft_sq = np.abs(np.fft.fft(x_normalized, axis=0))**2

    Y_fft_sq = Y_fft_sq[:nf, :, :]
    X_fft_sq = X_fft_sq[:nf, :, :]

    Svv = np.sum(Y_fft_sq, axis=(1, 2))
    Snm = np.sum(X_fft_sq, axis=(1, 2))

    denominator = Svv + Snm
    betaN = np.zeros_like(denominator, dtype=float)
    valid_den = denominator != 0
    betaN[valid_den] = Svv[valid_den] / denominator[valid_den]
    betaN[~valid_den] = np.nan

    return betaN

def lft_py(Y_fft_input, M_win, sig_bin_idx, N_channels):
    """
    Calcula a estatística para o Likelihood Function Test (LFT) - Detector baseado em potência.
    Y_fft_input: FFTs [Nfreqs, M_windows, N_channels]
    M_win: número de janelas (não usado diretamente no cálculo da estatística, mas para consistência)
    sig_bin_idx: índice do bin de frequência do sinal
    N_channels: número de canais
    """
    # A estatística LFT é a soma dos quadrados das magnitudes (potência) no bin do sinal
    # em todas as janelas e canais.
    lft_statistic = np.sum(np.abs(Y_fft_input[sig_bin_idx, :, :])**2)
    
    return lft_statistic

# --- Main Simulation Parameters ---
FS = 1000
NFFT_msc = 512
M_windows = 10
f_signal = 100 # Frequência do sinal de interesse
N_channels = 1 # Caso de 1 canal

target_pfa = 0.05
N_trials_H0 = 50000 # Para estimativa do limiar
N_trials_H1 = 5000  # Para estimativa de PD

snr_db_per_bin_range = np.arange(-30, 5 + 1, 1)

# --- Pre-calculations ---
total_samples_msc = NFFT_msc * M_windows
t_msc_segment = np.arange(NFFT_msc) / FS
hanning_window = np.hanning(NFFT_msc)

signal_bin_ideal_idx_msc = int(round(f_signal * NFFT_msc / FS))
f_signal_adjusted = signal_bin_ideal_idx_msc * FS / NFFT_msc
msc_signal_bin_py_idx = signal_bin_ideal_idx_msc

# --- Threshold Calculation (H0) for each method ---
thresholds = {}
msc_values_H0 = np.zeros(N_trials_H0)
mcsm_values_H0 = np.zeros(N_trials_H0)
lft_values_H0 = np.zeros(N_trials_H0)
mgbt_values_H0 = np.zeros(N_trials_H0)
sigma_n_sq_time = 1.0
print("Calculando limiares H0...")
for i in range(N_trials_H0):
    # Gerar ruído para todas as janelas e canais
    noise_all_windows_H0 = np.sqrt(sigma_n_sq_time) * np.random.randn(total_samples_msc, N_channels)
    
    # Preparar FFTs para MSC, MCSM, LFT
    # Y_fft_H0_all_windows será [Nfreqs, M_windows, N_channels]
    Y_fft_H0_all_windows = np.zeros((NFFT_msc // 2 + 1, M_windows, N_channels), dtype=complex)
    for k_win in range(M_windows):
        start_idx = k_win * NFFT_msc
        end_idx = (k_win + 1) * NFFT_msc
        for ch in range(N_channels):
            window_data = noise_all_windows_H0[start_idx:end_idx, ch]
            Y_fft_H0_all_windows[:, k_win, ch] = np.fft.rfft(window_data * hanning_window)
    
    # MSC Threshold (usando o primeiro canal)
    msc_spectrum_H0 = msc_fft_py(Y_fft_H0_all_windows[:, :, 0], M_windows)
    if msc_spectrum_H0 is not None and len(msc_spectrum_H0) > msc_signal_bin_py_idx:
        msc_values_H0[i] = msc_spectrum_H0[msc_signal_bin_py_idx]
    else: msc_values_H0[i] = np.nan

    # MCSM Threshold
    mcsm_spectrum_H0 = mcsm_py(noise_all_windows_H0, NFFT_msc, FS)
    if mcsm_spectrum_H0 is not None and len(mcsm_spectrum_H0) > msc_signal_bin_py_idx:
        mcsm_values_H0[i] = mcsm_spectrum_H0[msc_signal_bin_py_idx]
    else: mcsm_values_H0[i] = np.nan

    # LFT Threshold (soma dos quadrados das magnitudes no bin do sinal)
    lft_statistic_H0 = lft_py(Y_fft_H0_all_windows, M_windows, msc_signal_bin_py_idx, N_channels)
    lft_values_H0[i] = lft_statistic_H0

    # MGBT Threshold (y é ruído, x é outro ruído para comparação)
    noise_x_H0 = np.sqrt(sigma_n_sq_time) * np.random.randn(total_samples_msc, N_channels)
    mgbt_spectrum_H0 = mgbt_py(noise_all_windows_H0, noise_x_H0, NFFT_msc, FS)
    if mgbt_spectrum_H0 is not None and len(mgbt_spectrum_H0) > msc_signal_bin_py_idx:
        mgbt_values_H0[i] = mgbt_spectrum_H0[msc_signal_bin_py_idx]
    else: mgbt_values_H0[i] = np.nan

    if (i + 1) % (N_trials_H0 // 10) == 0:
        print(f'  H0 trial {i+1}/{N_trials_H0}')

# Calcula os percentis para os limiares
thresholds['MSC'] = scoreatpercentile(msc_values_H0[~np.isnan(msc_values_H0)], (1 - target_pfa) * 100)
thresholds['CSM'] = scoreatpercentile(mcsm_values_H0[~np.isnan(mcsm_values_H0)], (1 - target_pfa) * 100)
thresholds['LFT'] = scoreatpercentile(lft_values_H0[~np.isnan(lft_values_H0)], (1 - target_pfa) * 100)
thresholds['GBT'] = scoreatpercentile(mgbt_values_H0[~np.isnan(mgbt_values_H0)], (1 - target_pfa) * 100)

print("\nLimiares Calculados:")
for method, th in thresholds.items():
    print(f"  {method}: {th:.4f}")

# --- PD Calculation (H1) for each method ---
pd_results = {method: np.zeros(len(snr_db_per_bin_range)) for method in thresholds.keys()}

P_n_bin_expected = (2.0 / NFFT_msc) * sigma_n_sq_time

print("\nCalculando PD para diferentes SNRs (H1)...")
for snr_idx, current_snr_db in enumerate(snr_db_per_bin_range):
    current_snr_linear = 10**(current_snr_db / 10)
    P_s_bin_target = current_snr_linear * P_n_bin_expected
    A_signal = np.sqrt(2 * P_s_bin_target)
    
    detections = {method: 0 for method in thresholds.keys()}

    s_t_single_window = A_signal * np.sin(2 * np.pi * f_signal_adjusted * t_msc_segment)
    s_t_all_windows = np.tile(s_t_single_window, M_windows).reshape(-1, N_channels) # Reshape for N_channels
    s_t_all_windows = np.tile(s_t_single_window, M_windows).reshape(-1, N_channels)

    for i in range(N_trials_H1):
        # Gerar ruído para todas as janelas e canais
        n_t_all_windows = np.sqrt(sigma_n_sq_time) * np.random.randn(total_samples_msc, N_channels)
        x_t_all_windows = s_t_all_windows + n_t_all_windows
        
        # Preparar FFTs para todos os métodos
        Y_fft_H1_all_windows = np.zeros((NFFT_msc // 2 + 1, M_windows, N_channels), dtype=complex)
        for k_win in range(M_windows):
            start_idx = k_win * NFFT_msc
            end_idx = (k_win + 1) * NFFT_msc
            for ch in range(N_channels):
                window_data = x_t_all_windows[start_idx:end_idx, ch]
                Y_fft_H1_all_windows[:, k_win, ch] = np.fft.rfft(window_data * hanning_window)

        # MGBT precisa de um sinal de ruído separado 'x' para comparação
        n_t_all_windows_mgbt = np.sqrt(sigma_n_sq_time) * np.random.randn(total_samples_msc, N_channels)

        # Evaluate each method
        # Avaliar cada método
        msc_spectrum_H1 = msc_fft_py(Y_fft_H1_all_windows[:, :, 0], M_windows)
        if msc_spectrum_H1 is not None and len(msc_spectrum_H1) > msc_signal_bin_py_idx:
            if msc_spectrum_H1[msc_signal_bin_py_idx] > thresholds['MSC']:
                detections['MSC'] += 1

        mcsm_spectrum_H1 = mcsm_py(x_t_all_windows, NFFT_msc, FS)
        if mcsm_spectrum_H1 is not None and len(mcsm_spectrum_H1) > msc_signal_bin_py_idx:
            if mcsm_spectrum_H1[msc_signal_bin_py_idx] > thresholds['CSM']:
                detections['CSM'] += 1
        
        lft_statistic_H1 = np.sum(np.abs(Y_fft_H1_all_windows[msc_signal_bin_py_idx, :, 0])**2)
        lft_statistic_H1 = lft_py(Y_fft_H1_all_windows, M_windows, msc_signal_bin_py_idx, N_channels)
        if lft_statistic_H1 > thresholds['LFT']:
            detections['LFT'] += 1

        mgbt_spectrum_H1 = mgbt_py(x_t_all_windows, n_t_all_windows_mgbt, NFFT_msc, FS)
        if mgbt_spectrum_H1 is not None and len(mgbt_spectrum_H1) > msc_signal_bin_py_idx:
            if mgbt_spectrum_H1[msc_signal_bin_py_idx] > thresholds['GBT']:
                detections['GBT'] += 1
    
    # Calcular PD para cada método
    for method in thresholds.keys():
        pd_results[method][snr_idx] = detections[method] / N_trials_H1
    
    print(f"  SNR = {current_snr_db:.1f} dB: " + \
          ", ".join([f"PD({m}) = {pd_results[m][snr_idx]:.3f}" for m in thresholds.keys()]))

# --- Plotting ---
# --- Plotagem ---
plt.figure(figsize=(10, 6))
for method, pd_vals in pd_results.items():
    plt.plot(snr_db_per_bin_range, pd_vals, '-o', label=f'PD {method}')
print("\nSimulation complete.")
