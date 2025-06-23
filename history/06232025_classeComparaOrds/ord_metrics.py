import numpy as np
# No other imports like scipy.stats are strictly needed for these base metric calculations themselves,
# but they might be used if you were, for example, calculating theoretical distributions or p-values here.

# --- Univariate ORD Helper Functions ---

def msc_fft_py(Y_input_all_bins, M_win_actual):
    """
    Calculates Magnitude Squared Coherence for all bins.
    Args:
        Y_input_all_bins (np.ndarray): Complex FFT coefficients (num_bins, M_win_actual).
                                       If single bin, can be (1, M_win_actual).
        M_win_actual (int): Actual number of windows used.
    Returns:
        np.ndarray: MSC values for each bin (num_bins,).
    """
    if Y_input_all_bins.ndim == 1: # Handles case where a single bin's epochs are passed as 1D array
        Y_input_all_bins = Y_input_all_bins.reshape(1, -1)
        
    num_bins, M_win_data = Y_input_all_bins.shape

    if M_win_data == 0:
        return np.zeros(num_bins)
    if M_win_data < 2: # MSC is typically 0 or 1 (or undefined) for M_win < 2. Returning 0 for safety.
        return np.zeros(num_bins) 

    sum_Y_fft = np.sum(Y_input_all_bins, axis=1)
    numerator = np.abs(sum_Y_fft)**2

    sum_abs_Y_fft_sq = np.sum(np.abs(Y_input_all_bins)**2, axis=1)
    # The denominator uses M_win_data (actual windows) for scaling.
    denominator = M_win_data * sum_abs_Y_fft_sq

    msc = np.zeros(num_bins, dtype=float)
    # Avoid division by zero: where denominator is zero, msc remains zero.
    valid_den_mask = denominator != 0
    msc[valid_den_mask] = numerator[valid_den_mask] / denominator[valid_den_mask]
    
    msc[np.isnan(msc)] = 0.0 # Replace any NaNs that might arise (e.g., 0/0)
    msc[np.isinf(msc)] = 0.0 # Replace Infs (though less likely with this formulation)
    return msc

def csm_fft_py(Y_fft_all_windows, M_win_dummy=None): # M_win_dummy not used
    """
    Calculates Component Synchrony Measure (Phase Coherence) for all bins.
    Args:
        Y_fft_all_windows (np.ndarray): Complex FFT coefficients (num_bins, M_win_actual).
    Returns:
        np.ndarray: CSM values for each bin (num_bins,).
    """
    if Y_fft_all_windows.ndim == 1:
        Y_fft_all_windows = Y_fft_all_windows.reshape(1, -1)

    num_bins, M_win_actual = Y_fft_all_windows.shape
    if M_win_actual == 0:
        return np.zeros(num_bins)

    phases = np.angle(Y_fft_all_windows)
    mean_cos = np.mean(np.cos(phases), axis=1)
    mean_sin = np.mean(np.sin(phases), axis=1)
    csm_spectrum = mean_cos**2 + mean_sin**2
    csm_spectrum[np.isnan(csm_spectrum)] = 0.0
    return csm_spectrum

def lft_fft_py(Y_fft_all_windows_single_chan, M_win_dummy, L_neighbors):
    """
    Calculates Local F-Test for a single channel's full spectrum.
    Args:
        Y_fft_all_windows_single_chan (np.ndarray): Complex FFT (num_bins, M_win_actual).
        M_win_dummy (any): Not used, placeholder for consistent signature.
        L_neighbors (int): Number of neighboring bins on each side for noise estimation.
    Returns:
        np.ndarray: LFT values for each bin (num_bins,).
    """
    num_bins, M_win_actual = Y_fft_all_windows_single_chan.shape
    if M_win_actual == 0:
        return np.zeros(num_bins)
        
    power_spectrum_per_window = np.abs(Y_fft_all_windows_single_chan)**2
    avg_power_spectrum = np.mean(power_spectrum_per_window, axis=1)
    lft_spectrum = np.zeros(num_bins, dtype=float)

    for f_idx in range(num_bins):
        signal_power = avg_power_spectrum[f_idx]
        
        lower_indices = np.arange(max(0, f_idx - L_neighbors), f_idx)
        upper_indices = np.arange(f_idx + 1, min(num_bins, f_idx + L_neighbors + 1))
        all_neighbor_indices = np.unique(np.concatenate((lower_indices, upper_indices)).astype(int))
        valid_neighbors = all_neighbor_indices[(all_neighbor_indices >= 0) & (all_neighbor_indices < num_bins) & (all_neighbor_indices != f_idx)]
        
        if len(valid_neighbors) < L_neighbors / 2 : # Heuristic: need at least some distinct neighbors
             mean_noise_power = 1e-12 
        else:
            noise_power_bins = avg_power_spectrum[valid_neighbors]
            mean_noise_power = np.mean(noise_power_bins) if len(noise_power_bins) > 0 else 1e-12

        if mean_noise_power < 1e-12: lft_spectrum[f_idx] = signal_power / 1e-12
        else: lft_spectrum[f_idx] = signal_power / mean_noise_power
            
    lft_spectrum[np.isnan(lft_spectrum)] = 0.0
    # Cap Inf values based on the maximum finite LFT value or a large default
    finite_values = lft_spectrum[~np.isinf(lft_spectrum) & ~np.isnan(lft_spectrum)]
    cap_value = np.max(finite_values) if finite_values.size > 0 else 1e7 # Default large cap
    if cap_value == 0 and np.any(np.isinf(lft_spectrum)): cap_value = 1e7 # Ensure cap is large if all else is zero but infs exist
    lft_spectrum[np.isinf(lft_spectrum)] = cap_value
    return lft_spectrum

def gft_fft_py(Y_fft_all_windows_single_chan, target_bin_index, M_win_actual, noise_bins_mask=None, NFFT_val=512):
    """
    Calculates Global F-Test for a single channel.
    Args:
        Y_fft_all_windows_single_chan (np.ndarray): Complex FFT (num_bins, M_win_actual).
        target_bin_index (int): Index of the signal bin.
        M_win_actual (int): Actual number of windows used.
        noise_bins_mask (np.ndarray, optional): Boolean mask for noise bins.
        NFFT_val (int): FFT length, used if noise_bins_mask is None to exclude DC/Nyquist.
    Returns:
        float: GFT value.
    """
    num_bins, M_win_data = Y_fft_all_windows_single_chan.shape
    if M_win_data == 0: return np.nan
    
    power_spectrum_epochs = np.abs(Y_fft_all_windows_single_chan[:, :M_win_data])**2
    avg_power_spectrum = np.mean(power_spectrum_epochs, axis=1) 
    
    if not (0 <= target_bin_index < len(avg_power_spectrum)): return np.nan
    signal_power = avg_power_spectrum[target_bin_index]
    
    num_total_bins = len(avg_power_spectrum)
    if noise_bins_mask is None:
        noise_bins_mask = np.ones(num_total_bins, dtype=bool)
        noise_bins_mask[target_bin_index] = False
        if 0 < num_total_bins : noise_bins_mask[0] = False 
        if NFFT_val//2 < num_total_bins and NFFT_val//2 >= 0: # Ensure Nyquist index is valid
             noise_bins_mask[NFFT_val//2] = False 
        
    noise_power_values = avg_power_spectrum[noise_bins_mask]
    
    if len(noise_power_values) == 0: return np.nan 
    
    mean_noise_power = np.mean(noise_power_values)
    
    if mean_noise_power < 1e-12:
        return signal_power / 1e-12 
    gft_val = signal_power / mean_noise_power
    return gft_val if not np.isnan(gft_val) else 0.0

def ht2_fft_py(Y_fft_epoch_data_single_bin_M_win, M_win_actual):
    """
    Calculates Hotelling's T^2 for a single frequency bin across epochs.
    Args:
        Y_fft_epoch_data_single_bin_M_win (np.ndarray): Complex FFT coefficients for ONE bin, shape (M_win_actual,).
        M_win_actual (int): Number of windows/epochs.
    Returns:
        float: Hotelling's T^2 statistic.
    """
    if M_win_actual < 2 : return np.nan 
        
    Y_complex = Y_fft_epoch_data_single_bin_M_win.flatten()
    if len(Y_complex) != M_win_actual: return np.nan

    data_2d = np.array([Y_complex.real, Y_complex.imag]).T 
    
    try:
        mean_vector = np.mean(data_2d, axis=0) 
        if M_win_actual <= 2 : # For p=2, need M_win_actual > 2 for non-singular sample covariance
            # If M_win_actual is 2, cov matrix is often singular or T2 is not well defined in some stat packages.
            # Depending on the statistical definition, it might be 0 or a very large number.
            # For practical purposes in detection, if mean is non-zero it likely implies signal.
            # Let's return 0 if mean_vector is zero, else a large value, to indicate strong signal if M=2.
            return 1e7 if np.any(np.abs(mean_vector) > 1e-9) and M_win_actual==2 else 0.0
        
        # ddof=1 for sample covariance (N-1 in denominator)
        cov_matrix = np.cov(data_2d, rowvar=False, ddof=1) 
        
        if np.linalg.matrix_rank(cov_matrix) < 2:
             return 1e7 if np.any(np.abs(mean_vector) > 1e-9) else 0.0 

        inv_cov_matrix = np.linalg.inv(cov_matrix)
        t_squared = M_win_actual * (mean_vector.T @ inv_cov_matrix @ mean_vector)
        return t_squared if not np.isnan(t_squared) else 0.0
    except np.linalg.LinAlgError: 
        return 1e7 if np.any(np.abs(mean_vector) > 1e-9) else 0.0
    except Exception: 
        return np.nan

# --- "True" MORD Implementations ---

def mmsc_py(Y_fft_all_channels_all_windows_single_bin, M_win_actual):
    """
    Calculates Multiple Magnitude-Squared Coherence for a single frequency bin.
    Args:
        Y_fft_all_channels_all_windows_single_bin (np.ndarray): 
            Complex FFT coefficients. Shape: (num_channels, M_win_actual).
        M_win_actual (int): Actual number of windows used.
    Returns:
        float: MMSC value.
    """
    if Y_fft_all_channels_all_windows_single_bin.ndim == 1: # Passed single channel data
        # Reshape to (1, M_win_actual) and effectively calculate standard MSC
        Y_fft_all_channels_all_windows_single_bin = Y_fft_all_channels_all_windows_single_bin.reshape(1, -1)

    num_channels, M_win_data = Y_fft_all_channels_all_windows_single_bin.shape
    
    if M_win_data < 2: 
        return 0.0 

    sum_vectors_over_epochs = np.sum(Y_fft_all_channels_all_windows_single_bin, axis=1) # (num_channels,)
    numerator = np.linalg.norm(sum_vectors_over_epochs)**2 
                                                             
    sum_of_norm_sq_individual_epochs = 0.0
    for i in range(M_win_data):
        epoch_vector = Y_fft_all_channels_all_windows_single_bin[:, i] 
        sum_of_norm_sq_individual_epochs += np.linalg.norm(epoch_vector)**2
        
    if M_win_data == 0 or sum_of_norm_sq_individual_epochs < 1e-12: # Check for zero denominator
        return 0.0
        
    denominator = M_win_data * sum_of_norm_sq_individual_epochs
    
    mmsc_val = numerator / denominator
    return mmsc_val if not np.isnan(mmsc_val) else 0.0

def mcsm_py(Y_fft_all_channels_all_windows_single_bin, M_win_actual):
    """
    Calculates Multiple Component Synchrony Measure for a single frequency bin.
    Args:
        Y_fft_all_channels_all_windows_single_bin (np.ndarray): 
            Complex FFT. Shape: (num_channels, M_win_actual).
    Returns:
        float: MCSM value.
    """
    if Y_fft_all_channels_all_windows_single_bin.ndim == 1: # Single channel data
        Y_fft_all_channels_all_windows_single_bin = Y_fft_all_channels_all_windows_single_bin.reshape(1,-1)

    num_channels, M_win_data = Y_fft_all_channels_all_windows_single_bin.shape
    if M_win_data == 0: return 0.0

    phases_all_channels_epochs = np.angle(Y_fft_all_channels_all_windows_single_bin) 
    
    cos_phases_mc = np.cos(phases_all_channels_epochs) 
    sin_phases_mc = np.sin(phases_all_channels_epochs) 
    
    C_bar_epochs = np.mean(cos_phases_mc, axis=0) # Mean across channels for each epoch -> (M_win_data,)
    S_bar_epochs = np.mean(sin_phases_mc, axis=0) # -> (M_win_data,)
    
    mean_of_C_bar = np.mean(C_bar_epochs) # Mean across epochs
    mean_of_S_bar = np.mean(S_bar_epochs) 
    
    mcsm_val = mean_of_C_bar**2 + mean_of_S_bar**2
    return mcsm_val if not np.isnan(mcsm_val) else 0.0

def mlft_py(Y_fft_all_channels_all_windows_full_spectrum, target_bin_index, L_neighbors):
    """
    Calculates Multichannel Local F-Test for a target frequency bin.
    Args:
        Y_fft_all_channels_all_windows_full_spectrum (np.ndarray): 
            Complex FFT. Shape: (num_channels, num_freq_bins, M_win_actual).
        target_bin_index (int): Index of the signal bin.
        L_neighbors (int): Number of neighboring bins on each side.
    Returns:
        float: MLFT value for the target_bin_index.
    """
    if Y_fft_all_channels_all_windows_full_spectrum.ndim == 2: # Single channel full spectrum passed
        Y_fft_all_channels_all_windows_full_spectrum = Y_fft_all_channels_all_windows_full_spectrum.reshape(1, 
                                                                Y_fft_all_channels_all_windows_full_spectrum.shape[0], 
                                                                Y_fft_all_channels_all_windows_full_spectrum.shape[1])


    num_channels, num_freq_bins, M_win_actual = Y_fft_all_channels_all_windows_full_spectrum.shape
    if M_win_actual == 0: return np.nan

    power_mc_epochs_bins = np.abs(Y_fft_all_channels_all_windows_full_spectrum)**2 
    avg_power_across_channels_per_epoch = np.mean(power_mc_epochs_bins, axis=0) # (num_bins, M_win_actual)
    multichannel_avg_power_spectrum = np.mean(avg_power_across_channels_per_epoch, axis=1) # (num_bins,)

    if not (0 <= target_bin_index < num_freq_bins): return np.nan
    signal_power = multichannel_avg_power_spectrum[target_bin_index]
    
    lower_indices = np.arange(max(0, target_bin_index - L_neighbors), target_bin_index)
    upper_indices = np.arange(target_bin_index + 1, min(num_freq_bins, target_bin_index + L_neighbors + 1))
    all_neighbor_indices = np.unique(np.concatenate((lower_indices, upper_indices)).astype(int))
    valid_neighbors = all_neighbor_indices[(all_neighbor_indices >= 0) & (all_neighbor_indices < num_freq_bins) & (all_neighbor_indices != target_bin_index)]
    
    if len(valid_neighbors) < L_neighbors / 2: mean_noise_power = 1e-12
    else:
        noise_power_bins = multichannel_avg_power_spectrum[valid_neighbors]
        mean_noise_power = np.mean(noise_power_bins) if len(noise_power_bins) > 0 else 1e-12
    
    if mean_noise_power < 1e-12: mlft_val = signal_power / 1e-12
    else: mlft_val = signal_power / mean_noise_power
        
    mlft_val = 0.0 if np.isnan(mlft_val) else mlft_val
    if np.isinf(mlft_val):
        # Attempt to provide a reasonable cap if Inf occurs
        finite_spectrum_values = multichannel_avg_power_spectrum[~np.isinf(multichannel_avg_power_spectrum)]
        stable_denominator = mean_noise_power if mean_noise_power > 1e-12 else 1e-12
        if finite_spectrum_values.size > 0 and stable_denominator > 0 :
             mlft_val = np.max(finite_spectrum_values) / stable_denominator
        else: # Fallback if spectrum itself is problematic
             mlft_val = 1e7 # Default large value
    return mlft_val