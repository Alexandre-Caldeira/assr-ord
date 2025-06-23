import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import scoreatpercentile
import time
import pandas as pd
from collections import defaultdict
import pickle
import os

# --- Import from custom modules ---
from ord_metrics import (msc_fft_py, csm_fft_py, lft_fft_py, gft_fft_py, ht2_fft_py,
                       mmsc_py, mcsm_py, mlft_py)
from simulation_utils import (generate_correlated_noise_mc, wilson_score_interval_py,
                              adjust_freqs_get_bins_py)
from ord_processing import get_ord_value_enhanced, run_sequential_test_enhanced
from plotting_functions import plot_integrated_performance, plot_pd_vs_atl_sequential

# --- Simulation Parameters (Copied from your last main script example) ---
FS = 1000; NFFT = 512; L_neighbors_lft = 2
NUM_CHANNELS_SIM = 1
CHANNEL_CORRELATION = 0.0
if NUM_CHANNELS_SIM == 1: correlation_matrix_sim = np.array([[1.0]])
elif NUM_CHANNELS_SIM >= 2:
    if NUM_CHANNELS_SIM == 2: correlation_matrix_sim = np.array([[1.0, CHANNEL_CORRELATION], [CHANNEL_CORRELATION, 1.0]])
    else:
        correlation_matrix_sim = np.eye(NUM_CHANNELS_SIM)
        for i in range(NUM_CHANNELS_SIM):
            for j in range(i + 1, NUM_CHANNELS_SIM):
                correlation_matrix_sim[i,j] = correlation_matrix_sim[j,i] = CHANNEL_CORRELATION**(abs(i-j))
print(f"Simulating with {NUM_CHANNELS_SIM} channel(s), Correlation Matrix:\n{correlation_matrix_sim}")
signal_frequencies_nominal = np.array([82])
fp_reference_frequencies_nominal = np.array([41, 51, 61, 71])
noise_only_frequencies_H0_secondary = np.array([101, 111, 121, 131])
noise_only_frequencies_H1_evaluation = np.array([30, 50, 70, 120])
target_pfa = 0.05
N_trials_H0 = 200  # REDUCED FOR TESTING
N_trials_H1 = 50   # REDUCED FOR TESTING
M_MAX_WINDOWS_OPTIONS = [20, 60] # Example, use [10, 30, 60, 90] for full run
snr_db_per_bin_range = np.arange(-30, 5 + 1, 10) # Coarser for faster test
NDC_OPTIONS = [1, 15, 30] # Example, use list(range(1,31)) for full run
hanning_win = np.hanning(NFFT)
sigma_n_sq_time = 1.0
P_n_bin_expected = (2.0 / NFFT) * sigma_n_sq_time
RESULTS_DIR = "simulation_results_full"
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

# --- ORD Configurations (Copied from your last main script example) ---
ord_configurations_to_run = []
for base_ord_name in ['MSC', 'CSM', 'LFT', 'GFT', 'HT2']:
    config = {'id': f'{base_ord_name}_ch1_f0', 'base_ord': base_ord_name, 'num_channels_for_ord': 1, 'q_num_harmonics': 0}
    if base_ord_name == 'LFT': config['L_lft'] = L_neighbors_lft
    ord_configurations_to_run.append(config)
# Add more configurations as needed based on your previous setup...

# --- Main Simulation ---
all_results_list = []
signal_frequencies_adj, signal_bin_indices_py = adjust_freqs_get_bins_py(signal_frequencies_nominal, NFFT, FS)
fp_reference_frequencies_adj, fp_ref_bin_indices_py = adjust_freqs_get_bins_py(fp_reference_frequencies_nominal, NFFT, FS)
noise_only_H0_sec_adj, noise_only_H0_sec_bins = adjust_freqs_get_bins_py(noise_only_frequencies_H0_secondary, NFFT, FS)
noise_only_H1_eval_adj, noise_only_H1_eval_bins = adjust_freqs_get_bins_py(noise_only_frequencies_H1_evaluation, NFFT, FS)

target_signal_bin_fundamental = signal_bin_indices_py[0] if len(signal_bin_indices_py) > 0 else -1
if target_signal_bin_fundamental == -1: raise ValueError("No signal frequency defined!")

start_total_sim_time = time.time()

for M_MAX_CURRENT in M_MAX_WINDOWS_OPTIONS:
    print(f"\n===== M_MAX_WINDOWS = {M_MAX_CURRENT} =====")
    
    for ord_config in ord_configurations_to_run:
        current_ord_id = ord_config['id']
        print(f"\n--- ORD Config: {current_ord_id} (Base: {ord_config.get('base_ord','N/A')}) ---")

        # --- H0 Phase ---
        print(f"    H0 phase for {current_ord_id} (M_MAX={M_MAX_CURRENT})...")
        h0_ord_values_at_fp_ref_dist_all_trials = [] 
        h0_max_ord_across_fp_ref_trials = np.full(N_trials_H0, np.nan)
        h0_ord_values_at_noise_secondary_dist_all_trials = []

        for i_h0 in range(N_trials_H0):
            noise_mc_h0 = generate_correlated_noise_mc(NFFT * M_MAX_CURRENT, NUM_CHANNELS_SIM, sigma_n_sq_time, correlation_matrix_sim)
            Y_fft_H0_trial_all_ch = np.zeros((NUM_CHANNELS_SIM, NFFT // 2 + 1, M_MAX_CURRENT), dtype=complex)
            for ch_idx in range(NUM_CHANNELS_SIM):
                for k_win in range(M_MAX_CURRENT):
                    start_fft = k_win * NFFT; end_fft = (k_win + 1) * NFFT
                    Y_fft_H0_trial_all_ch[ch_idx, :, k_win] = np.fft.rfft(noise_mc_h0[start_fft:end_fft, ch_idx] * hanning_win)
            
            current_trial_ords_at_fp_ref = []
            for fp_fund_bin in fp_ref_bin_indices_py:
                ord_val = get_ord_value_enhanced(Y_fft_H0_trial_all_ch, M_MAX_CURRENT, ord_config, fp_fund_bin, FS, NFFT)
                if not np.isnan(ord_val): current_trial_ords_at_fp_ref.append(ord_val)
            if current_trial_ords_at_fp_ref:
                h0_ord_values_at_fp_ref_dist_all_trials.extend(current_trial_ords_at_fp_ref) # Store all individual values
                h0_max_ord_across_fp_ref_trials[i_h0] = np.max(current_trial_ords_at_fp_ref)
            
            current_trial_ords_at_noise_secondary = []
            for ns_bin in noise_only_H0_sec_bins:
                ord_val_ns = get_ord_value_enhanced(Y_fft_H0_trial_all_ch, M_MAX_CURRENT, ord_config, ns_bin, FS, NFFT)
                if not np.isnan(ord_val_ns): current_trial_ords_at_noise_secondary.append(ord_val_ns)
            if current_trial_ords_at_noise_secondary: 
                h0_ord_values_at_noise_secondary_dist_all_trials.extend(current_trial_ords_at_noise_secondary)

        valid_max_ord_H0 = h0_max_ord_across_fp_ref_trials[~np.isnan(h0_max_ord_across_fp_ref_trials)]
        if len(valid_max_ord_H0) < N_trials_H0 * 0.1:
            print(f"    Warning: Insufficient valid H0 trials for threshold for {current_ord_id} (found {len(valid_max_ord_H0)}). Skipping this ORD config for M_MAX={M_MAX_CURRENT}.")
            continue
        ord_threshold_current_config = scoreatpercentile(valid_max_ord_H0, (1 - target_pfa) * 100)
        
        pfa_h0_primary_val, pfa_h0_primary_ci = wilson_score_interval_py(np.sum(valid_max_ord_H0 > ord_threshold_current_config), len(valid_max_ord_H0))
        tnr_h0_primary_val = 1.0 - pfa_h0_primary_val if not np.isnan(pfa_h0_primary_val) else np.nan
        print(f"      Threshold (Single Test): {ord_threshold_current_config:.4f} (Emp. PFA_H0_primary: {pfa_h0_primary_val:.4f})")

        h0_secondary_detections = np.sum(np.array(h0_ord_values_at_noise_secondary_dist_all_trials) > ord_threshold_current_config)
        h0_secondary_total_evals = len(h0_ord_values_at_noise_secondary_dist_all_trials)
        pfa_h0_secondary_val, pfa_h0_secondary_ci = wilson_score_interval_py(h0_secondary_detections, h0_secondary_total_evals)
        tnr_h0_secondary_val = 1.0 - pfa_h0_secondary_val if not np.isnan(pfa_h0_secondary_val) else np.nan
        print(f"      Emp. PFA_H0_secondary (on other noise freqs): {pfa_h0_secondary_val:.4f}")

        result_H0_single = {
            'ORD_ID': current_ord_id, 'M_MAX': M_MAX_CURRENT, 'SNR_dB': 'H0_Overall', 'NDC': 'Single', 
            'Test_Type': 'Single', 'Threshold': ord_threshold_current_config,
            'PFA_H0_primary': pfa_h0_primary_val, 'PFA_H0_primary_CI_low':pfa_h0_primary_ci[0], 'PFA_H0_primary_CI_high':pfa_h0_primary_ci[1],
            'TNR_H0_primary': tnr_h0_primary_val,
            'PFA_H0_secondary': pfa_h0_secondary_val, 'PFA_H0_secondary_CI_low':pfa_h0_secondary_ci[0], 'PFA_H0_secondary_CI_high':pfa_h0_secondary_ci[1],
            'TNR_H0_secondary': tnr_h0_secondary_val,
            'H0_fp_ref_dist': list(h0_ord_values_at_fp_ref_dist_all_trials),
            'H0_noise_secondary_dist': list(h0_ord_values_at_noise_secondary_dist_all_trials)
        }
        for k_nan in ['PD_TPR', 'FNR', 'PFA_H1_noise', 'TNR_H1_noise', 'ATL_ANW']: result_H0_single[k_nan] = np.nan
        all_results_list.append(result_H0_single)
        
        print(f"    Re-evaluating H0 for sequential tests with threshold={ord_threshold_current_config:.4f}...")
        h0_seq_pfa_results = defaultdict(lambda: {'detections': 0, 'atl_list': []})
        # This part requires generating H0 data again, or storing it if memory allows. For N_trials_H0, regenerating is safer.
        for i_h0_seq in range(N_trials_H0):
            noise_mc_h0_seq = generate_correlated_noise_mc(NFFT * M_MAX_CURRENT, NUM_CHANNELS_SIM, sigma_n_sq_time, correlation_matrix_sim)
            Y_fft_H0_trial_all_ch_seq = np.zeros((NUM_CHANNELS_SIM, NFFT // 2 + 1, M_MAX_CURRENT), dtype=complex)
            for ch_idx in range(NUM_CHANNELS_SIM):
                for k_win in range(M_MAX_CURRENT):
                    start_fft = k_win * NFFT; end_fft = (k_win + 1) * NFFT
                    Y_fft_H0_trial_all_ch_seq[ch_idx, :, k_win] = np.fft.rfft(noise_mc_h0_seq[start_fft:end_fft, ch_idx] * hanning_win)
            
            target_h0_seq_fund_bin = fp_ref_bin_indices_py[0] # Simplified: use first fp_ref_bin for seq H0 PFA
            for ndc_val in NDC_OPTIONS:
                m_min_curr = max(1, ndc_val); base_ord = ord_config.get('base_ord')
                if base_ord in ['MSC', 'HT2', 'MMSC', 'MCSM'] and m_min_curr < 2: m_min_curr = 2
                if NDC_OPTIONS and ndc_val > m_min_curr : m_min_curr = ndc_val


                detected_seq, windows_used = run_sequential_test_enhanced(
                    Y_fft_H0_trial_all_ch_seq, ord_config, ord_threshold_current_config, 
                    target_h0_seq_fund_bin, m_min_curr, M_MAX_CURRENT, ndc_val, FS, NFFT
                )
                if detected_seq: h0_seq_pfa_results[ndc_val]['detections'] += 1
                h0_seq_pfa_results[ndc_val]['atl_list'].append(windows_used)
        
        for ndc_val in NDC_OPTIONS:
            pfa_s_val, pfa_s_ci = wilson_score_interval_py(h0_seq_pfa_results[ndc_val]['detections'], N_trials_H0)
            atl_s_h0_val = np.mean(h0_seq_pfa_results[ndc_val]['atl_list']) if h0_seq_pfa_results[ndc_val]['atl_list'] else np.nan
            print(f"      NDC={ndc_val}: PFA_seq(H0)={pfa_s_val:.4f}, Avg_ATL_seq(H0)={atl_s_h0_val:.2f}")
            result_H0_seq_ndc = {
                'ORD_ID': current_ord_id, 'M_MAX': M_MAX_CURRENT, 'SNR_dB': 'H0_Overall_Seq', 'NDC': ndc_val,
                'Test_Type': 'Sequential', 'Threshold': ord_threshold_current_config,
                'PFA_H0_primary': pfa_s_val, 'PFA_H0_primary_CI_low':pfa_s_ci[0], 'PFA_H0_primary_CI_high':pfa_s_ci[1],
                'TNR_H0_primary': 1.0-pfa_s_val if not np.isnan(pfa_s_val) else np.nan,
                'ATL_ANW': atl_s_h0_val,
            }
            for k_nan_seq in ['PFA_H0_secondary','TNR_H0_secondary','PD_TPR', 'FNR', 'PFA_H1_noise', 'TNR_H1_noise', 'H0_fp_ref_dist', 'H0_noise_secondary_dist', 'H1_signal_dist', 'H1_noise_eval_dist']: result_H0_seq_ndc[k_nan_seq] = np.nan
            all_results_list.append(result_H0_seq_ndc)

        # --- H1 Phase ---
        print(f"    H1 phase for {current_ord_id} (M_MAX={M_MAX_CURRENT})...")
        for snr_db in snr_db_per_bin_range:
            current_snr_linear = 10**(snr_db / 10)
            P_s_bin_target = current_snr_linear * P_n_bin_expected
            A_signal = np.sqrt(2 * P_s_bin_target) if P_s_bin_target >= 0 else 0 # Ensure A_signal is non-negative

            h1_detections_single_count = 0
            # Store raw ORD values from H1 trials at signal frequency
            h1_ord_values_at_signal_dist_snr = [] 
            # Store raw ORD values from H1 trials at H1_evaluation noise bins
            h1_ord_values_at_noise_eval_dist_snr = []
            
            h1_seq_pd_results = defaultdict(lambda: {'detections': 0, 'atl_list': []})
            h1_seq_pfa_noise_bins_results = defaultdict(lambda: {'fp_detections':0, 'total_noise_bin_trials':0})


            for i_h1 in range(N_trials_H1):
                s_t_mc = np.zeros((NFFT * M_MAX_CURRENT, NUM_CHANNELS_SIM))
                if A_signal > 0 and target_signal_bin_fundamental != -1:
                    t_total_h1 = np.arange(NFFT * M_MAX_CURRENT) / FS
                    f0_hz_h1 = target_signal_bin_fundamental * FS / NFFT
                    base_signal_f0_h1 = A_signal * np.sin(2 * np.pi * f0_hz_h1 * t_total_h1)
                    for ch_idx in range(NUM_CHANNELS_SIM): s_t_mc[:, ch_idx] = base_signal_f0_h1
                
                noise_mc_h1 = generate_correlated_noise_mc(NFFT * M_MAX_CURRENT, NUM_CHANNELS_SIM, sigma_n_sq_time, correlation_matrix_sim)
                x_t_mc_h1 = s_t_mc + noise_mc_h1
                
                Y_fft_H1_trial_all_ch = np.zeros((NUM_CHANNELS_SIM, NFFT // 2 + 1, M_MAX_CURRENT), dtype=complex)
                for ch_idx in range(NUM_CHANNELS_SIM):
                    for k_win in range(M_MAX_CURRENT):
                        start_fft = k_win * NFFT; end_fft = (k_win + 1) * NFFT
                        Y_fft_H1_trial_all_ch[ch_idx, :, k_win] = np.fft.rfft(x_t_mc_h1[start_fft:end_fft, ch_idx] * hanning_win)

                # Single Test H1 - PD
                ord_val_h1_single = get_ord_value_enhanced(Y_fft_H1_trial_all_ch, M_MAX_CURRENT, ord_config, 
                                                           target_signal_bin_fundamental, FS, NFFT)
                if not np.isnan(ord_val_h1_single):
                    h1_ord_values_at_signal_dist_snr.append(ord_val_h1_single)
                    if ord_val_h1_single > ord_threshold_current_config:
                        h1_detections_single_count += 1
                
                # Single Test H1 - PFA on noise_only_frequencies_H1_evaluation
                for h1_noise_bin in noise_only_H1_eval_bins:
                    ord_val_h1_noise = get_ord_value_enhanced(Y_fft_H1_trial_all_ch, M_MAX_CURRENT, ord_config,
                                                              h1_noise_bin, FS, NFFT) 
                    if not np.isnan(ord_val_h1_noise): 
                        h1_ord_values_at_noise_eval_dist_snr.append(ord_val_h1_noise)
                
                # Sequential Tests H1
                for ndc_val in NDC_OPTIONS:
                    m_min_curr = max(1, ndc_val); base_ord = ord_config.get('base_ord')
                    if base_ord in ['MSC', 'HT2', 'MMSC', 'MCSM'] and m_min_curr < 2: m_min_curr = 2
                    if NDC_OPTIONS and ndc_val > m_min_curr : m_min_curr = ndc_val

                    detected_seq, windows_used = run_sequential_test_enhanced(
                        Y_fft_H1_trial_all_ch, ord_config, ord_threshold_current_config, 
                        target_signal_bin_fundamental, m_min_curr, M_MAX_CURRENT, ndc_val, FS, NFFT
                    )
                    if detected_seq: h1_seq_pd_results[ndc_val]['detections'] += 1
                    h1_seq_pd_results[ndc_val]['atl_list'].append(windows_used)
                    
                    # PFA_H1_noise Sequential (simplified: one noise bin)
                    target_h1_noise_seq_bin = noise_only_H1_eval_bins[0] if len(noise_only_H1_eval_bins)>0 else -1
                    if target_h1_noise_seq_bin != -1:
                        detected_seq_noise, _ = run_sequential_test_enhanced(
                            Y_fft_H1_trial_all_ch, ord_config, ord_threshold_current_config,
                            target_h1_noise_seq_bin, m_min_curr, M_MAX_CURRENT, ndc_val, FS, NFFT)
                        if detected_seq_noise: h1_seq_pfa_noise_bins_results[ndc_val]['fp_detections'] += 1
                        h1_seq_pfa_noise_bins_results[ndc_val]['total_noise_bin_trials'] +=1
            
            pd_single_val, pd_single_ci = wilson_score_interval_py(h1_detections_single_count, N_trials_H1)
            fnr_single_val = 1.0 - pd_single_val if not np.isnan(pd_single_val) else np.nan
            
            h1_noise_single_fp_count = np.sum(np.array(h1_ord_values_at_noise_eval_dist_snr) > ord_threshold_current_config)
            pfa_h1_noise_single_val, pfa_h1_noise_single_ci = wilson_score_interval_py(h1_noise_single_fp_count, len(h1_ord_values_at_noise_eval_dist_snr) if h1_ord_values_at_noise_eval_dist_snr else 0)
            tnr_h1_noise_single_val = 1.0 - pfa_h1_noise_single_val if not np.isnan(pfa_h1_noise_single_val) else np.nan

            all_results_list.append({
                'ORD_ID': current_ord_id, 'M_MAX': M_MAX_CURRENT, 'SNR_dB': snr_db, 'NDC': 'Single', 
                'Test_Type': 'Single', 'Threshold': ord_threshold_current_config,
                'PD_TPR': pd_single_val, 'PD_CI_low':pd_single_ci[0],'PD_CI_high':pd_single_ci[1], 'FNR': fnr_single_val,
                'PFA_H1_noise': pfa_h1_noise_single_val, 'PFA_H1_noise_CI_low':pfa_h1_noise_single_ci[0],'PFA_H1_noise_CI_high':pfa_h1_noise_single_ci[1],
                'TNR_H1_noise': tnr_h1_noise_single_val, 'ATL_ANW': M_MAX_CURRENT,
                'H1_signal_dist': list(h1_ord_values_at_signal_dist_snr), 
                'H1_noise_eval_dist': list(h1_ord_values_at_noise_eval_dist_snr),
                'PFA_H0_primary':np.nan, 'TNR_H0_primary':np.nan,'PFA_H0_secondary':np.nan,'TNR_H0_secondary':np.nan,
                 'H0_fp_ref_dist': [], 'H0_noise_secondary_dist': []
            })

            for ndc_val in NDC_OPTIONS:
                pd_s_val, pd_s_ci = wilson_score_interval_py(h1_seq_pd_results[ndc_val]['detections'], N_trials_H1)
                fnr_s_val = 1.0 - pd_s_val if not np.isnan(pd_s_val) else np.nan
                atl_s_h1_val = np.mean(h1_seq_pd_results[ndc_val]['atl_list']) if h1_seq_pd_results[ndc_val]['atl_list'] else np.nan
                
                pfa_h1_noise_s_val, pfa_h1_noise_s_ci = wilson_score_interval_py(
                    h1_seq_pfa_noise_bins_results[ndc_val]['fp_detections'], 
                    h1_seq_pfa_noise_bins_results[ndc_val]['total_noise_bin_trials']
                )
                tnr_h1_noise_s_val = 1.0 - pfa_h1_noise_s_val if not np.isnan(pfa_h1_noise_s_val) else np.nan

                all_results_list.append({
                    'ORD_ID': current_ord_id, 'M_MAX': M_MAX_CURRENT, 'SNR_dB': snr_db, 'NDC': ndc_val, 
                    'Test_Type': 'Sequential', 'Threshold': ord_threshold_current_config,
                    'PD_TPR': pd_s_val, 'PD_CI_low':pd_s_ci[0],'PD_CI_high':pd_s_ci[1], 'FNR': fnr_s_val,
                    'PFA_H1_noise': pfa_h1_noise_s_val, 'PFA_H1_noise_CI_low':pfa_h1_noise_s_ci[0],'PFA_H1_noise_CI_high':pfa_h1_noise_s_ci[1],
                    'TNR_H1_noise': tnr_h1_noise_s_val, 'ATL_ANW': atl_s_h1_val,
                    'PFA_H0_primary':np.nan, 'TNR_H0_primary':np.nan,'PFA_H0_secondary':np.nan,'TNR_H0_secondary':np.nan,
                    'H0_fp_ref_dist': [], 'H0_noise_secondary_dist': [], 'H1_signal_dist': [], 'H1_noise_eval_dist': []
                })
            print(f"      SNR={snr_db:.1f}dB processed: PD_single={pd_single_val:.3f}, PD_seq(NDC={NDC_OPTIONS[0]})={h1_seq_pd_results[NDC_OPTIONS[0]]['detections']/N_trials_H1 if N_trials_H1 > 0 else np.nan:.3f}, Avg_ATL(NDC={NDC_OPTIONS[0]})={np.mean(h1_seq_pd_results[NDC_OPTIONS[0]]['atl_list']) if h1_seq_pd_results[NDC_OPTIONS[0]]['atl_list'] else np.nan:.2f}")

# --- Finalize and Save Results ---
# [The rest of the script: DataFrame creation, saving, plotting calls as before]
# ... (Copied from previous response)
# --- Save Full Results ---
results_df = pd.DataFrame(all_results_list) if all_results_list else pd.DataFrame()
if not results_df.empty:
    print("\n--- Simulation Results DataFrame Head (Example with Dummy Data) ---")
    print(results_df.head())
    results_filename = os.path.join(RESULTS_DIR, f"example_main_sim_MMAX{M_MAX_WINDOWS_OPTIONS[0]}_NCh{NUM_CHANNELS_SIM}_Corr{CHANNEL_CORRELATION}.pkl")
    try:
        with open(results_filename, 'wb') as f:
            pickle.dump(results_df, f)
        print(f"Full results saved to {results_filename}")
    except Exception as e:
        print(f"Error saving results: {e}")
else:
    print("Warning: No results were generated to save or plot.")

end_total_sim_time = time.time()
print(f"\nTotal Simulation Time (Example with Dummy H0/H1): {(end_total_sim_time - start_total_sim_time)/60:.2f} minutes")


# --- Save Full Results ---
results_df = pd.DataFrame(all_results_list) if all_results_list else pd.DataFrame()
if not results_df.empty:
    print("\n--- Simulation Results DataFrame Head ---")
    print(results_df.head())
    results_filename = os.path.join(RESULTS_DIR, f"full_main_sim_MMAX{M_MAX_WINDOWS_OPTIONS[0]}_NCh{NUM_CHANNELS_SIM}_Corr{CHANNEL_CORRELATION}.pkl")
    try:
        with open(results_filename, 'wb') as f:
            pickle.dump(results_df, f)
        print(f"Full results saved to {results_filename}")
    except Exception as e:
        print(f"Error saving results: {e}")
else:
    print("Warning: No results were generated to save or plot.")

end_total_sim_time = time.time()
print(f"\nTotal Simulation Time: {(end_total_sim_time - start_total_sim_time)/60:.2f} minutes")

# --- Plotting Calls (using functions from plotting_functions.py) ---
if not results_df.empty:
    plotable_ord_ids = [oid for oid in results_df['ORD_ID'].unique() if isinstance(oid, str) and 'H0_Overall' not in oid] # Ensure string and not H0 summary
    if plotable_ord_ids:
        first_ord_id_to_plot = plotable_ord_ids[0]
        m_max_values_for_ord = results_df[results_df['ORD_ID'] == first_ord_id_to_plot]['M_MAX'].unique()
        if m_max_values_for_ord.size > 0:
            first_m_max_to_plot = m_max_values_for_ord[0]
            print(f"\n--- Generating Integrated Plot for: {first_ord_id_to_plot}, M_MAX={first_m_max_to_plot} ---")
            plot_integrated_performance(results_df, first_ord_id_to_plot, first_m_max_to_plot, target_pfa, NDC_OPTIONS)
            print(f"\n--- Generating PD vs ATL Plot for: {first_ord_id_to_plot}, M_MAX={first_m_max_to_plot} ---")
            plot_pd_vs_atl_sequential(results_df, first_ord_id_to_plot, first_m_max_to_plot, target_pfa, NDC_OPTIONS, M_MAX_WINDOWS_OPTIONS)
        else: print(f"No M_MAX data found for ORD_ID: {first_ord_id_to_plot} to generate plots.")
    else: print("No plottable ORD_IDs found in results (after filtering H0 summary rows).")
else: print("No data in results_df to plot.")

print("\nMain Simulation Script with Full H0/H1 Logic - Execution Complete.")