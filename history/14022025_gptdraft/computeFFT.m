function fftSignals = computeFFT(signals, params) 
% computeFFT computes the FFT along each window (row-wise) and 
% returns only the positive frequencies. 
    fftSignals = fft(signals); 
    fftSignals = fftSignals(:, 1:params.nBins,:); 
end
