import numpy as np
from scipy.stats import norm # For wilson_score_interval_py

def generate_correlated_noise_mc(num_samples_total, num_channels, sigma_sq_time, correlation_matrix):
    """
    Generates multi-channel correlated Gaussian noise.
    Args:
        num_samples_total (int): Total number of time samples to generate for each channel.
        num_channels (int): Number of channels.
        sigma_sq_time (float): Variance of noise for each channel (assumed equal before correlation).
        correlation_matrix (np.ndarray): Target (num_channels x num_channels) correlation matrix.
                                         Must be positive semi-definite.
    Returns:
        np.ndarray: Correlated noise array of shape (num_samples_total, num_channels).
    """
    if num_channels == 0:
        return np.array([]).reshape(num_samples_total, 0)
    if correlation_matrix.shape[0] != num_channels or correlation_matrix.shape[1] != num_channels:
        raise ValueError("Correlation matrix dimensions must match num_channels.")

    std_dev = np.sqrt(sigma_sq_time)
    # Create a diagonal matrix of standard deviations
    D = np.diag(np.full(num_channels, std_dev))
    # Covariance matrix C = D @ R @ D, where R is the correlation matrix
    covariance_matrix = D @ correlation_matrix @ D
    
    # Generate independent standard normal noise (num_channels x num_samples_total for matrix multiplication)
    independent_noise = np.random.randn(num_channels, num_samples_total)
    
    try:
        # Cholesky decomposition: Cov = L @ L.T
        # L is a lower triangular matrix
        L = np.linalg.cholesky(covariance_matrix)
        # Correlated noise Y = L @ X (where X is standard normal)
        # Result is (num_channels x num_samples_total)
        correlated_noise_transposed = L @ independent_noise
        # Transpose to get (num_samples_total x num_channels)
        return correlated_noise_transposed.T
    except np.linalg.LinAlgError:
        # This can happen if covariance_matrix is not positive definite
        # print("Warning: Cholesky decomposition failed for noise generation. Using independent noise as fallback.")
        # Fallback to independent noise if Cholesky fails
        return np.sqrt(sigma_sq_time) * np.random.randn(num_samples_total, num_channels)

def wilson_score_interval_py(k, n, confidence_level=0.95):
    """
    Calculates the Wilson score interval for a binomial proportion.
    Args:
        k (int or float): Number of successes.
        n (int or float): Total number of trials.
        confidence_level (float): The desired confidence level (e.g., 0.95 for 95% CI).
    Returns:
        tuple: (phat, (lower_bound, upper_bound))
               phat is the observed proportion (k/n).
               The second element is a tuple for the confidence interval.
               Returns (np.nan, (np.nan, np.nan)) if input is invalid (e.g., n=0).
    """
    if n == 0 or np.isnan(n) or np.isnan(k) or k < 0 or n < k:
        return np.nan, (np.nan, np.nan)
    
    k = float(np.nan_to_num(k)) # Ensure k is a float, handle potential NaNs passed
    n = float(n)                # Ensure n is a float

    alpha = 1 - confidence_level
    z = norm.ppf(1 - alpha / 2) # z-score for the confidence level

    phat = k / n
    
    # Wilson score interval formula components
    term1_numerator = phat + (z**2) / (2 * n)
    
    sqrt_arg = (phat * (1 - phat) / n) + (z**2) / (4 * n**2)
    # Ensure sqrt_arg is non-negative due to potential floating point inaccuracies
    if sqrt_arg < 0: 
        sqrt_arg = 0 
    term2_numerator = z * np.sqrt(sqrt_arg)
    
    denominator = 1 + (z**2) / n
    
    lower_bound = (term1_numerator - term2_numerator) / denominator
    upper_bound = (term1_numerator + term2_numerator) / denominator
    
    # Ensure bounds are within [0, 1]
    pci_l = max(0.0, lower_bound)
    pci_u = min(1.0, upper_bound)
    
    return phat, (pci_l, pci_u)

def adjust_freqs_get_bins_py(nominal_freqs, NFFT_val, FS_rate):
    """
    Adjusts nominal frequencies to the center of FFT bins and returns bin indices.
    Args:
        nominal_freqs (np.ndarray or list): Array of nominal frequencies.
        NFFT_val (int): FFT length.
        FS_rate (int or float): Sampling frequency.
    Returns:
        tuple: (adjusted_frequencies (np.ndarray), bin_indices (np.ndarray of int))
               Bin indices are for use with np.fft.rfft output (0 to NFFT//2).
    """
    if not isinstance(nominal_freqs, np.ndarray):
        nominal_freqs = np.array(nominal_freqs)
        
    adj_freqs = np.zeros_like(nominal_freqs, dtype=float)
    bin_indices = np.zeros_like(nominal_freqs, dtype=int)
    
    freq_resolution = FS_rate / NFFT_val
    
    for i, f_val in enumerate(nominal_freqs):
        # Calculate the ideal (possibly fractional) bin index
        ideal_bin_idx_float = f_val / freq_resolution
        # Round to the nearest integer bin index
        bin_idx_int = int(round(ideal_bin_idx_float))
        
        # Ensure bin index is within the valid range for rfft output (0 to NFFT//2)
        bin_idx_int = max(0, min(bin_idx_int, NFFT_val // 2))
        
        adj_freqs[i] = bin_idx_int * freq_resolution
        bin_indices[i] = bin_idx_int
        
    return adj_freqs, bin_indices

# You could add other general utilities here later, for example:
# - Functions for specific signal generation patterns (e.g., AM tones, chirps)
# - More advanced noise generation (e.g., colored noise)
# - Data loading/saving helpers if they become complex