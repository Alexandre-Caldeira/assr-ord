function MSCvalues = computeMSC(fftSignals, params) 
% computeMSC calculates a simplified magnitude-squared coherence (MSC) value per stage. 
    K = params.K;
    totalWindows = params.duration; 
    windowsPerStage = floor(totalWindows / K); 

    % MSCvalues = zeros(1, K); 
    MSCvalues = zeros(params.nChannels,K,params.nBins);

    for k = 1:K 
        idxStart = (k-1)*windowsPerStage + 1; 
        idxEnd = idxStart + windowsPerStage - 1; 
        
        for channel = 1:params.nChannels
            MSCvalues(channel,k,:) = msc_fft( ...
                squeeze(fftSignals(channel, :, idxStart:idxEnd)), ...
                windowsPerStage);

        end
    end 
end
