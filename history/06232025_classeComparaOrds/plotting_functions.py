import numpy as np
import matplotlib.pyplot as plt
import pandas as pd # Used for handling the DataFrame input to plots

def plot_integrated_performance(df_results, ord_id_to_plot, m_max_to_plot, target_pfa_val, ndc_options_for_h0_lines):
    """
    Generates an "Integrated Performance Plot" for a specific ORD configuration.
    Shows PD/TPR, FNR vs. SNR, and various PFA/TNR lines based on H0 data.

    Args:
        df_results (pd.DataFrame): The main DataFrame containing all simulation results.
        ord_id_to_plot (str): The unique ID of the ORD configuration to plot.
        m_max_to_plot (int): The M_MAX value for which to plot.
        target_pfa_val (float): The target PFA used in simulations.
        ndc_options_for_h0_lines (list): List of NDC values, used to find a fallback H0 PFA line
                                         if single test H0 data is missing.
    """
    
    # Filter for H0 overall summary (single test) for this ORD_ID and M_MAX
    h0_summary_df_single = df_results[
        (df_results['ORD_ID'] == ord_id_to_plot) &
        (df_results['M_MAX'] == m_max_to_plot) &
        (df_results['SNR_dB'] == 'H0_Overall') &  # Specific key for single test H0 summary
        (df_results['Test_Type'] == 'Single')
    ]
    
    # Initialize PFA/TNR values to NaN or fallbacks
    pfa_h0_pri_line = target_pfa_val # Fallback to target PFA
    tnr_h0_pri_line = 1.0 - target_pfa_val
    pfa_h0_sec_line = np.nan
    tnr_h0_sec_line = np.nan

    if not h0_summary_df_single.empty:
        pfa_h0_pri_line = h0_summary_df_single['PFA_H0_primary'].iloc[0]
        tnr_h0_pri_line = h0_summary_df_single['TNR_H0_primary'].iloc[0]
        pfa_h0_sec_line = h0_summary_df_single['PFA_H0_secondary'].iloc[0]
        tnr_h0_sec_line = h0_summary_df_single['TNR_H0_secondary'].iloc[0]
    else:
        # Fallback: If no 'H0_Overall' for Single, try 'H0_Overall_Seq' for the first NDC
        # This is a heuristic if the H0 summary for single test wasn't stored under that exact key.
        first_ndc = ndc_options_for_h0_lines[0] if ndc_options_for_h0_lines else 1
        h0_summary_df_seq_fallback = df_results[
            (df_results['ORD_ID'] == ord_id_to_plot) &
            (df_results['M_MAX'] == m_max_to_plot) &
            (df_results['SNR_dB'] == 'H0_Overall_Seq') & # Key for sequential H0 summary
            (df_results['NDC'] == first_ndc)
        ]
        if not h0_summary_df_seq_fallback.empty:
            pfa_h0_pri_line = h0_summary_df_seq_fallback['PFA_H0_primary'].iloc[0]
            tnr_h0_pri_line = h0_summary_df_seq_fallback['TNR_H0_primary'].iloc[0]
            # PFA_H0_secondary might not be available for seq H0 summary, so it remains NaN or from target
        # else:
            # print(f"Warning: No H0 summary data found for PFA lines for {ord_id_to_plot}, M_MAX={m_max_to_plot}.")


    # Data for PD/FN/PFA_H1 vs SNR (single test for this ORD_ID and M_MAX)
    h1_plot_df_single = df_results[
        (df_results['ORD_ID'] == ord_id_to_plot) &
        (df_results['M_MAX'] == m_max_to_plot) &
        (pd.to_numeric(df_results['SNR_dB'], errors='coerce').notna()) & # Ensures SNR_dB is numeric (H1 data)
        (df_results['Test_Type'] == 'Single')
    ].copy() # Use .copy() to avoid SettingWithCopyWarning

    if not h1_plot_df_single.empty:
        h1_plot_df_single['SNR_dB_num'] = pd.to_numeric(h1_plot_df_single['SNR_dB'])
        h1_plot_df_single = h1_plot_df_single.sort_values(by='SNR_dB_num')
    else:
        print(f"No H1 data (Single Test) to plot for {ord_id_to_plot} at M_MAX {m_max_to_plot}.")
        # Optionally, still create an empty plot or return
        # For now, we'll proceed and it will plot only H0 lines if no H1 data.

    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot PD (TPR) and FNR from H1 single test data
    if not h1_plot_df_single.empty:
        ax.plot(h1_plot_df_single['SNR_dB_num'], h1_plot_df_single['PD_TPR'], marker='o', linestyle='-', 
                label=f'PD (TPR) - {ord_id_to_plot}')
        if 'PD_CI_low' in h1_plot_df_single.columns and 'PD_CI_high' in h1_plot_df_single.columns:
             ax.fill_between(h1_plot_df_single['SNR_dB_num'], 
                             h1_plot_df_single['PD_CI_low'].fillna(0), # fillna for cases where CI might be nan
                             h1_plot_df_single['PD_CI_high'].fillna(1), alpha=0.2)
        ax.plot(h1_plot_df_single['SNR_dB_num'], h1_plot_df_single['FNR'], marker='x', linestyle='--', 
                label=f'FNR - {ord_id_to_plot}')

    # Plot PFA/TNR lines
    ax.axhline(pfa_h0_pri_line, color='red', linestyle=':', 
               label=f'PFA_H0_primary ({pfa_h0_pri_line:.3f})')
    if not np.isnan(tnr_h0_pri_line):
        ax.axhline(tnr_h0_pri_line, color='blue', linestyle=':', 
                   label=f'TNR_H0_primary ({tnr_h0_pri_line:.3f})')
    
    if not np.isnan(pfa_h0_sec_line):
        ax.axhline(pfa_h0_sec_line, color='salmon', linestyle='-.', 
                   label=f'PFA_H0_secondary ({pfa_h0_sec_line:.3f})')
    if not np.isnan(tnr_h0_sec_line):
        ax.axhline(tnr_h0_sec_line, color='skyblue', linestyle='-.', 
                   label=f'TNR_H0_secondary ({tnr_h0_sec_line:.3f})')

    # Average PFA_H1_noise and TNR_H1_noise
    if not h1_plot_df_single.empty:
        avg_pfa_h1_noise = h1_plot_df_single['PFA_H1_noise'].mean()
        avg_tnr_h1_noise = h1_plot_df_single['TNR_H1_noise'].mean()
        if not np.isnan(avg_pfa_h1_noise):
            ax.axhline(avg_pfa_h1_noise, color='lightcoral', linestyle='--', 
                       label=f'Avg PFA_H1_noise ({avg_pfa_h1_noise:.3f})')
        if not np.isnan(avg_tnr_h1_noise):
            ax.axhline(avg_tnr_h1_noise, color='lightblue', linestyle='--', 
                       label=f'Avg TNR_H1_noise ({avg_tnr_h1_noise:.3f})')

    ax.set_title(f'Integrated Performance (Single Test): {ord_id_to_plot}\n(M_MAX={m_max_to_plot}, Target PFA ≈ {target_pfa_val:.2f})')
    ax.set_xlabel('SNR per Bin (dB)')
    ax.set_ylabel('Rate / Probability')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
    ax.grid(True)
    ax.set_ylim([-0.05, 1.05])
    plt.tight_layout(rect=[0,0,0.78,1]) # Adjust for legend outside
    plt.show()


def plot_pd_vs_atl_sequential(df_results, ord_id_to_plot, m_max_to_plot, target_pfa_val, ndc_options_plot_arg, M_MAX_WINDOWS_OPTIONS_ARG_PLOT):
    """
    Plots PD vs. SNR and ATL vs. SNR for sequential tests, comparing different NDCs
    and the single test performance for a specific ORD configuration.
    """
    fig, axs = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    fig.suptitle(f"Sequential vs. Single Test: {ord_id_to_plot}, M_MAX={m_max_to_plot}\n(Target PFA ≈ {target_pfa_val:.2f})", fontsize=14)

    # --- PD Plot (axs[0]) ---
    ax1 = axs[0]
    # Single test H1 data
    single_h1_df = df_results[
        (df_results['ORD_ID'] == ord_id_to_plot) &
        (df_results['M_MAX'] == m_max_to_plot) &
        (pd.to_numeric(df_results['SNR_dB'], errors='coerce').notna()) &
        (df_results['Test_Type'] == 'Single')
    ].copy()

    if not single_h1_df.empty:
        single_h1_df['SNR_dB_num'] = pd.to_numeric(single_h1_df['SNR_dB'])
        single_h1_df = single_h1_df.sort_values(by='SNR_dB_num')
        ax1.plot(single_h1_df['SNR_dB_num'], single_h1_df['PD_TPR'], marker='x', linestyle='--', 
                 label='PD Single Test', color='black')

    # Sequential test H1 data for various NDCs
    for ndc_val in ndc_options_plot_arg:
        seq_h1_df_ndc = df_results[
            (df_results['ORD_ID'] == ord_id_to_plot) &
            (df_results['M_MAX'] == m_max_to_plot) &
            (pd.to_numeric(df_results['SNR_dB'], errors='coerce').notna()) &
            (df_results['Test_Type'] == 'Sequential') & 
            (df_results['NDC'] == ndc_val) # Ensure NDC is numeric for comparison
        ].copy()
        if not seq_h1_df_ndc.empty:
            seq_h1_df_ndc['SNR_dB_num'] = pd.to_numeric(seq_h1_df_ndc['SNR_dB'])
            seq_h1_df_ndc = seq_h1_df_ndc.sort_values(by='SNR_dB_num')
            ax1.plot(seq_h1_df_ndc['SNR_dB_num'], seq_h1_df_ndc['PD_TPR'], marker='o', 
                     label=f'PD Seq (NDC={ndc_val})')
    
    ax1.set_ylabel('Probability of Detection (PD/TPR)')
    ax1.set_title('Detection Performance (PD vs SNR)')
    ax1.legend(loc='lower right', fontsize='small'); ax1.grid(True); ax1.set_ylim([-0.05, 1.05])

    # --- ATL Plot (axs[1]) ---
    ax2 = axs[1]
    if not single_h1_df.empty: # ATL for single test under H1 is always M_MAX
         ax2.plot(single_h1_df['SNR_dB_num'], single_h1_df['ATL_ANW'], marker='x', linestyle='--', 
                  label='ATL Single Test', color='black')

    for ndc_val in ndc_options_plot_arg:
        seq_h1_df_ndc = df_results[
            (df_results['ORD_ID'] == ord_id_to_plot) &
            (df_results['M_MAX'] == m_max_to_plot) &
            (pd.to_numeric(df_results['SNR_dB'], errors='coerce').notna()) &
            (df_results['Test_Type'] == 'Sequential') & 
            (df_results['NDC'] == ndc_val)
        ].copy()
        if not seq_h1_df_ndc.empty:
            seq_h1_df_ndc['SNR_dB_num'] = pd.to_numeric(seq_h1_df_ndc['SNR_dB'])
            seq_h1_df_ndc = seq_h1_df_ndc.sort_values(by='SNR_dB_num')
            ax2.plot(seq_h1_df_ndc['SNR_dB_num'], seq_h1_df_ndc['ATL_ANW'], marker='o', 
                     label=f'ATL Seq (NDC={ndc_val})')
    
    ax2.set_xlabel('SNR per Bin (dB)')
    ax2.set_ylabel('Average Test Length (Windows)')
    ax2.set_title('Test Efficiency (ATL vs SNR)')
    ax2.legend(loc='upper right', fontsize='small'); ax2.grid(True);
    # Dynamic Ylim based on M_MAX or actual max ATL
    max_atl_val = m_max_to_plot
    if M_MAX_WINDOWS_OPTIONS_ARG_PLOT and len(M_MAX_WINDOWS_OPTIONS_ARG_PLOT)>0:
         max_atl_val = M_MAX_WINDOWS_OPTIONS_ARG_PLOT[0] # Assuming it's a list/array

    df_atl_max_check = df_results[ (df_results['ORD_ID'] == ord_id_to_plot) & (df_results['M_MAX'] == m_max_to_plot) & (pd.to_numeric(df_results['SNR_dB'], errors='coerce').notna())]['ATL_ANW']
    if not df_atl_max_check.empty and pd.notnull(df_atl_max_check.max()):
        plot_upper_ylim = max(max_atl_val, df_atl_max_check.max()) * 1.1
    else:
        plot_upper_ylim = max_atl_val * 1.1

    ax2.set_ylim([0, plot_upper_ylim if plot_upper_ylim > 0 else m_max_to_plot + 5])
    
    plt.tight_layout(rect=[0, 0, 1, 0.95]) # Adjust for suptitle
    plt.show()