import numpy as np

# Import the specific ORD functions from ord_metrics.py
# Ensure ord_metrics.py is in the same directory or accessible via PYTHONPATH
from ord_metrics import (
    msc_fft_py, csm_fft_py, lft_fft_py, gft_fft_py, ht2_fft_py,
    mmsc_py, mcsm_py, mlft_py
)

def get_ord_value_enhanced(Y_fft_trial_data_all_channels,  # Shape: (num_sim_channels, num_freq_bins, current_m_windows)
                           current_m_windows,
                           ord_config,  # Dictionary defining the ORD
                           target_fundamental_bin_index,
                           FS_rate, NFFT_val):
    """
    Calculates the final ORD statistic based on the ord_config.
    This function handles:
    - Selecting the base ORD (univariate or true MORD).
    - Applying q-sample logic (combining results from fundamental + harmonics).
    - Applying channel combination for derived MORDs (e.g., aMSC, pCSM).
    """
    if current_m_windows == 0:
        return np.nan

    base_ord_name = ord_config['base_ord']
    num_sim_channels_available = Y_fft_trial_data_all_channels.shape[0]
    # num_channels_for_ord specifies how many of the available channels this ORD config should use.
    num_ch_to_use_for_ord = min(ord_config.get('num_channels_for_ord', num_sim_channels_available), num_sim_channels_available)

    # Data to be processed by this specific ORD configuration
    # Shape: (num_ch_to_use_for_ord, num_freq_bins, current_m_windows)
    Y_data_this_ord = Y_fft_trial_data_all_channels[:num_ch_to_use_for_ord, :, :current_m_windows]

    # Determine all frequency bins to evaluate (fundamental + harmonics)
    harmonic_bins_to_evaluate = [target_fundamental_bin_index]
    if ord_config.get('q_num_harmonics', 0) > 0:
        fundamental_freq_hz = target_fundamental_bin_index * FS_rate / NFFT_val
        for h_idx in range(1, ord_config['q_num_harmonics'] + 1):
            harmonic_freq_ideal = fundamental_freq_hz * (h_idx + 1) # f0, 2*f0, 3*f0...
            harmonic_bin = int(round(harmonic_freq_ideal * NFFT_val / FS_rate))
            # Ensure harmonic bin is valid and not already included
            if 0 <= harmonic_bin <= NFFT_val // 2 and harmonic_bin not in harmonic_bins_to_evaluate:
                harmonic_bins_to_evaluate.append(harmonic_bin)
    harmonic_bins_to_evaluate = sorted(list(set(harmonic_bins_to_evaluate)))

    # --- Logic branches based on ORD type ---

    # Case 1: True MORDs (MMSC, MCSM, MLFT) - these inherently combine channel info.
    # q-sample logic is applied to their outputs if specified.
    if base_ord_name in ['MMSC', 'MCSM', 'MLFT']:
        ord_stats_for_harmonics = []
        for bin_idx in harmonic_bins_to_evaluate:
            if not (0 <= bin_idx < Y_data_this_ord.shape[1]): continue # Skip invalid bin index

            ord_val_one_harmonic = np.nan
            if base_ord_name == 'MMSC':
                # mmsc_py expects (num_ch_to_use, M_win) for a single bin
                ord_val_one_harmonic = mmsc_py(Y_data_this_ord[:, bin_idx, :], current_m_windows)
            elif base_ord_name == 'MCSM':
                # mcsm_py expects (num_ch_to_use, M_win) for a single bin
                ord_val_one_harmonic = mcsm_py(Y_data_this_ord[:, bin_idx, :], current_m_windows)
            elif base_ord_name == 'MLFT':
                # mlft_py expects (num_ch_to_use, num_bins, M_win) and the target bin
                ord_val_one_harmonic = mlft_py(Y_data_this_ord, bin_idx, ord_config.get('L_lft', 2))
            
            if not np.isnan(ord_val_one_harmonic):
                ord_stats_for_harmonics.append(ord_val_one_harmonic)

        if not ord_stats_for_harmonics: return np.nan
        # Apply q-combination to the MORD outputs from different harmonics
        if ord_config.get('q_num_harmonics', 0) == 0 or len(ord_stats_for_harmonics) == 1:
            return ord_stats_for_harmonics[0]
        else:
            q_method = ord_config.get('q_combination', 'sum')
            valid_stats = [s for s in ord_stats_for_harmonics if not np.isnan(s)]
            if not valid_stats: return np.nan
            if q_method == 'sum': return np.sum(valid_stats)
            elif q_method == 'average': return np.mean(valid_stats)
            elif q_method == 'product': return np.prod([max(s, 1e-9) for s in valid_stats]) # Avoid log(0) if using geometric mean idea
            elif q_method == 'max': return np.max(valid_stats)
            else: return np.nan

    # Case 2: Univariate base ORDs (MSC, CSM, LFT, GFT, HT2)
    # These are calculated per channel. If num_ch_to_use_for_ord > 1 AND channel_combination is specified,
    # we calculate the (potentially q-sampled) ORD for each channel, then combine these channel values.
    else:
        per_channel_final_ord_values = [] # Stores the (q-combined) ORD value for each channel
        for chan_idx in range(num_ch_to_use_for_ord):
            Y_single_chan_data = Y_data_this_ord[chan_idx, :, :] # (num_freq_bins, current_m_windows)
            
            ord_stats_for_harmonics_this_channel = []
            for bin_idx in harmonic_bins_to_evaluate:
                if not (0 <= bin_idx < Y_single_chan_data.shape[0]): continue

                ord_stat_one_bin = np.nan
                # Data for single bin, single channel, all epochs for MSC, CSM, HT2: (current_m_windows,)
                # Data for full spectrum, single channel, all epochs for LFT, GFT: (num_freq_bins, current_m_windows)
                
                if base_ord_name == 'MSC':
                    ord_stat_one_bin = msc_fft_py(Y_single_chan_data[bin_idx, :].reshape(1, current_m_windows), current_m_windows)[0]
                elif base_ord_name == 'CSM':
                    ord_stat_one_bin = csm_fft_py(Y_single_chan_data[bin_idx, :].reshape(1, current_m_windows), current_m_windows)[0]
                elif base_ord_name == 'HT2':
                    ord_stat_one_bin = ht2_fft_py(Y_single_chan_data[bin_idx, :], current_m_windows)
                elif base_ord_name == 'LFT':
                    lft_spectrum_this_chan = lft_fft_py(Y_single_chan_data, current_m_windows, ord_config.get('L_lft', 2))
                    if bin_idx < len(lft_spectrum_this_chan): ord_stat_one_bin = lft_spectrum_this_chan[bin_idx]
                elif base_ord_name == 'GFT':
                    ord_stat_one_bin = gft_fft_py(Y_single_chan_data, bin_idx, current_m_windows, NFFT_val=NFFT_val)
                
                if not np.isnan(ord_stat_one_bin):
                    ord_stats_for_harmonics_this_channel.append(ord_stat_one_bin)
            
            # q-Combination for this channel
            if not ord_stats_for_harmonics_this_channel:
                per_channel_final_ord_values.append(np.nan)
                continue
            
            q_val_this_channel = np.nan
            if ord_config.get('q_num_harmonics', 0) == 0 or len(ord_stats_for_harmonics_this_channel) == 1:
                q_val_this_channel = ord_stats_for_harmonics_this_channel[0]
            else:
                q_method = ord_config.get('q_combination', 'sum')
                valid_stats = [s for s in ord_stats_for_harmonics_this_channel if not np.isnan(s)]
                if not valid_stats: q_val_this_channel = np.nan
                elif q_method == 'sum': q_val_this_channel = np.sum(valid_stats)
                elif q_method == 'average': q_val_this_channel = np.mean(valid_stats)
                elif q_method == 'product': q_val_this_channel = np.prod([max(s, 1e-9) for s in valid_stats])
                elif q_method == 'max': q_val_this_channel = np.max(valid_stats)
            per_channel_final_ord_values.append(q_val_this_channel)

        # Now, apply channel combination if specified and relevant
        valid_per_channel_ords = [s for s in per_channel_final_ord_values if not np.isnan(s)]
        if not valid_per_channel_ords: return np.nan

        if num_ch_to_use_for_ord == 1 or not ord_config.get('channel_combination'):
            return valid_per_channel_ords[0] # Result from the single channel used
        else: # Derived MORD - apply channel combination
            combo_method = ord_config['channel_combination']
            if combo_method == 'average': return np.mean(valid_per_channel_ords)
            elif combo_method == 'product': # Geometric mean
                prods = [max(s, 1e-9) for s in valid_per_channel_ords]
                return np.prod(prods) ** (1.0 / len(prods)) if prods else np.nan
            elif combo_method == 'max': return np.max(valid_per_channel_ords)
            else: return np.nan # Unknown combination

    return np.nan # Should not be reached if logic is complete


def run_sequential_test_enhanced(Y_fft_trial_data_all_channels, ord_config_dict, threshold,
                                 target_fundamental_bin, M_MIN_param, M_MAX, NDC, FS_rate, NFFT_val_seq):
    """
    Runs a single sequential test trial using the enhanced ORD calculation.
    Args:
        Y_fft_trial_data_all_channels (np.ndarray): (num_sim_channels, num_freq_bins, M_MAX_for_trial_data)
        ord_config_dict (dict): Configuration for the ORD to be used.
        threshold (float): Detection threshold for the ORD.
        target_fundamental_bin (int): Bin index of the fundamental frequency to test.
        M_MIN_param (int): Minimum number of windows before first decision.
        M_MAX (int): Maximum number of windows for the trial.
        NDC (int): Number of consecutive detections required.
        FS_rate (int): Sampling frequency.
        NFFT_val_seq (int): FFT length used for processing.
    Returns:
        tuple: (detected (bool), num_windows_used (int))
    """
    consecutive_detections = 0
    
    # Determine actual minimum windows based on ORD type and NDC
    actual_m_min = M_MIN_param
    base_ord = ord_config_dict.get('base_ord')
    if base_ord in ['MSC', 'HT2', 'MMSC', 'MCSM'] and actual_m_min < 2: # These need at least 2 epochs
        actual_m_min = 2
    actual_m_min = max(actual_m_min, 1) # Ensure M_MIN is at least 1
    if NDC > actual_m_min : # Can't have NDC detections before M_MIN if M_MIN is larger
        actual_m_min = NDC # Effectively, M_MIN for decision making is max(configured_M_MIN, NDC, ord_specific_min)
    
    for m_current in range(1, M_MAX + 1):
        if m_current < actual_m_min:
            # If m_current is less than NDC, we can't satisfy the NDC condition yet, so reset counter
            if m_current < NDC: 
                consecutive_detections = 0
            continue # Skip decision if below effective M_MIN for decision making

        # Pass data up to m_current for processing
        current_epoch_data = Y_fft_trial_data_all_channels[:, :, :m_current]
        
        ord_val = get_ord_value_enhanced(current_epoch_data, 
                                         m_current, # This is the number of windows being used now
                                         ord_config_dict,
                                         target_fundamental_bin, 
                                         FS_rate, NFFT_val_seq)
        
        if not np.isnan(ord_val) and ord_val > threshold:
            consecutive_detections += 1
        else:
            consecutive_detections = 0 # Reset if below threshold or NaN
        
        if consecutive_detections >= NDC:
            return True, m_current # Detected
            
    return False, M_MAX # Not detected within M_MAX windows