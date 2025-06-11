classdef ORDTester < handle
    % Applies statistical tests (e.g., Beta CGST) to MSC data.

    properties (SetAccess = protected)
        desiredAlpha = 0.05;
        thresholdCache
        signalFrequencies
        noiseFrequencies
        allTestFrequencies
        nSignalFrequencies
        nNoiseFrequencies
        nAllTestFrequencies
        noiseFlagIndex
    end

    methods
        function obj = ORDTester(signalFrequencies, noiseFrequencies, options)
            arguments
                signalFrequencies (1,:) {mustBeNumeric, mustBePositive}
                noiseFrequencies (1,:) {mustBeNumeric, mustBePositive} % Allow empty
                options.desiredAlpha (1,1) {mustBeNumeric, ...
                mustBeGreaterThan(options.desiredAlpha,0), ...
                mustBeLessThan(options.desiredAlpha,1)} = 0.05
            end

            if numel(signalFrequencies) ~= numel(unique(signalFrequencies)), error('ORDTester:InvalidInput', 'Input signalFrequencies must contain unique values.'); end
            if ~isempty(noiseFrequencies) && numel(noiseFrequencies) ~= numel(unique(noiseFrequencies)), error('ORDTester:InvalidInput', 'Input noiseFrequencies must contain unique values.'); end

            obj.signalFrequencies = sort(signalFrequencies);
            obj.noiseFrequencies = sort(noiseFrequencies);
            obj.desiredAlpha = options.desiredAlpha;
            obj.allTestFrequencies = sort([obj.signalFrequencies, obj.noiseFrequencies]);

            commonFreqs = intersect(obj.signalFrequencies, obj.noiseFrequencies);
            if ~isempty(commonFreqs), error('Signal frequencies and noise frequencies must be distinct. Common: %s', mat2str(commonFreqs)); end

            obj.nSignalFrequencies = numel(obj.signalFrequencies);
            obj.nNoiseFrequencies = numel(obj.noiseFrequencies);
            obj.nAllTestFrequencies = numel(obj.allTestFrequencies);

            if obj.nNoiseFrequencies > 0
                [member, location] = ismember(min(obj.noiseFrequencies), obj.allTestFrequencies);
                 if ~member, error('Could not determine noise frequency start index.'); end
                 obj.noiseFlagIndex = location;
            else
                obj.noiseFlagIndex = obj.nAllTestFrequencies + 1;
            end

            obj.thresholdCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
        end

        function [decisionResults, stageAlphas, stageGammas] = runBetaCGST(obj, mscData, M, K)
             arguments
                obj
                mscData (:,:,:) {mustBeNumeric}
                M (1,1) {mustBeInteger, mustBeGreaterThanOrEqual(M, 2)} % Corrected validation
                K (1,1) {mustBeInteger, mustBePositive}
            end

            if K ~= size(mscData, 2), error('K (%d) does not match MSC data dim 2 (%d).', K, size(mscData, 2)); end

            [stageAlphas, stageGammas] = obj.getBetaCGSTThresholds(M, K); % Pass correct M

            freqIndices = obj.allTestFrequencies + 1;
            maxBinIndex = size(mscData, 1);
            if any(freqIndices > maxBinIndex)
                warning('Clipping test frequencies (%s) exceeding max bin index (%d).', mat2str(obj.allTestFrequencies(freqIndices > maxBinIndex)), maxBinIndex);
                freqIndices = min(freqIndices, maxBinIndex);
                freqIndices = max(freqIndices, 1);
            end
            mscAtTestFreqs = mscData(freqIndices, :, :);

            nChannels = size(mscData, 3);
            TP = zeros(obj.nAllTestFrequencies, K, nChannels);
            TN = zeros(obj.nAllTestFrequencies, K, nChannels);
            FP = zeros(obj.nAllTestFrequencies, K, nChannels);
            FN = zeros(obj.nAllTestFrequencies, K, nChannels);

            cumulativeMSC = cumsum(mscAtTestFreqs, 2);

            isSignal = false(obj.nAllTestFrequencies, 1);
            if obj.noiseFlagIndex > 1, isSignal(1 : obj.noiseFlagIndex - 1) = true; end
            isNoise = ~isSignal;
            if obj.noiseFlagIndex <= obj.nAllTestFrequencies, isNoise(obj.noiseFlagIndex : end) = true; isSignal(obj.noiseFlagIndex : end) = false; end

            for ch = 1:nChannels
                decided_ch = false(obj.nAllTestFrequencies, 1);
                for k_stage = 1:K
                    currentCumulativeMSC_ch = cumulativeMSC(:, k_stage, ch);
                    detected = currentCumulativeMSC_ch > stageAlphas(k_stage);
                    stopped = currentCumulativeMSC_ch <= stageGammas(k_stage);
                    eligible = ~decided_ch;

                    TP(eligible & isSignal & detected, k_stage, ch) = 1;
                    FN(eligible & isSignal & stopped, k_stage, ch) = 1;
                    FP(eligible & isNoise & detected, k_stage, ch) = 1;
                    TN(eligible & isNoise & stopped, k_stage, ch) = 1;

                    decided_now = eligible & (detected | stopped);
                    decided_ch(decided_now) = true;

                     if k_stage == K
                         final_undecided_ch = ~decided_ch;
                         if any(final_undecided_ch)
                             FN(final_undecided_ch & isSignal, k_stage, ch) = 1;
                             TN(final_undecided_ch & isNoise, k_stage, ch) = 1;
                             decided_ch(final_undecided_ch) = true;
                         end
                     end
                end
            end

            decisionResults = struct('TP', TP, 'TN', TN, 'FP', FP, 'FN', FN, ...
                                     'allTestFrequencies', obj.allTestFrequencies, ...
                                     'signalFrequencies', obj.signalFrequencies, ...
                                     'noiseFrequencies', obj.noiseFrequencies);
        end


        function [stageAlphas, stageGammas] = getBetaCGSTThresholds(obj, M, K)
             arguments
                obj
                M (1,1) {mustBeInteger, mustBeGreaterThanOrEqual(M, 2)} % Corrected validation
                K (1,1) {mustBeInteger, mustBePositive}
            end
            cacheKey = sprintf('M%dK%d_A%.4f', M, K, obj.desiredAlpha);

            if obj.thresholdCache.isKey(cacheKey)
                cachedData = obj.thresholdCache(cacheKey);
                stageAlphas = cachedData.stageAlphas; stageGammas = cachedData.stageGammas;
            else
                alpha = obj.desiredAlpha;
                Alpha_k = ones(1, K) * (alpha / K);
                Gamma_k = ((1 - alpha) / K) .* ones(1, K);

                Resolution = 1e6; % Using higher resolution with FFT opt.
                Xvalues = 0 : (1/Resolution) : 1;
                targetLen = Resolution + 1;

                 if M == 2, Null = ones(1, targetLen); else, Null = betapdf(Xvalues, 1, M - 1); end
                 if sum(Null) < eps, error('Initial Beta PDF calculation failed for M=%d', M); end
                 Null = Null / sum(Null);
                 if length(Null) ~= targetLen, error('Initial Null PDF length mismatch'); end
                Chi2_Norm = Null;

                aThresholds = zeros(1, K); gThresholds = zeros(1, K);
                aThresholds(1) = 1 - Alpha_k(1).^(1 ./ (M - 1));
                gThresholds(1) = 1 - (1 - Gamma_k(1)).^(1 ./ (M - 1));
                aThresholds(1) = max(0, min(1, aThresholds(1)));
                gThresholds(1) = max(0, min(1, gThresholds(1)));

                currentNull = Null;
                probMassAccounted = Alpha_k(1) + Gamma_k(1);

                fftLen = 2 * targetLen - 1;
                fft_Chi2_Norm = fft(Chi2_Norm, fftLen);

                for k = 2:K
                    TruncInd_Ra_prev = round(aThresholds(k-1) * Resolution) + 1;
                    TruncInd_Rg_prev = round(gThresholds(k-1) * Resolution) + 1;
                    NullTrunc = currentNull;
                    if TruncInd_Ra_prev <= targetLen, NullTrunc(TruncInd_Ra_prev:end) = 0; end
                    if TruncInd_Rg_prev >= 1, NullTrunc(1:TruncInd_Rg_prev) = 0; end

                    sumTrunc = sum(NullTrunc);
                    if sumTrunc < 1e-12, warning('Null distribution truncated to zero at stage k=%d.', k); aThresholds(k:end) = 1; gThresholds(k:end) = 0; break; end
                    NullTrunc = NullTrunc / sumTrunc;

                    fft_NullTrunc = fft(NullTrunc, fftLen);
                    conv_result_fft = ifft(fft_NullTrunc .* fft_Chi2_Norm, fftLen);
                    Null2_unscaled = real(conv_result_fft(1:targetLen));
                    Null2_unscaled(Null2_unscaled < 0) = 0;

                    remainingProbMass = 1 - probMassAccounted;
                    currentSum = sum(Null2_unscaled);
                    if currentSum < 1e-12 || remainingProbMass < 1e-12, warning('Normalization factor issue (sum=%.2e, rem=%.2e) at stage k=%d.',currentSum, remainingProbMass, k); aThresholds(k:end) = 1; gThresholds(k:end) = 0; break; end
                    Null2 = Null2_unscaled * (remainingProbMass / currentSum);

                    targetAlphaMass = Alpha_k(k); cumulativeSumRight = cumsum(Null2(end:-1:1)); idxAlpha = find(cumulativeSumRight >= targetAlphaMass, 1, 'first');
                    if isempty(idxAlpha), aThresholds(k) = 1; else, aThresholds(k) = (targetLen - idxAlpha) / Resolution; end

                    targetGammaMass = Gamma_k(k); cumulativeSumLeft = cumsum(Null2); idxGamma = find(cumulativeSumLeft >= targetGammaMass, 1, 'first');
                    if isempty(idxGamma), gThresholds(k) = 0; else, gThresholds(k) = (idxGamma - 1) / Resolution; end

                    aThresholds(k) = min(max(aThresholds(k), 0), 1); gThresholds(k) = min(max(gThresholds(k), 0), 1);
                    if aThresholds(k) <= gThresholds(k), gThresholds(k) = max(0, aThresholds(k) - 1/Resolution); end
                    if aThresholds(k) < gThresholds(k), warning('Threshold calc error: Alpha < Gamma at k=%d.', k); gThresholds(k) = aThresholds(k); end

                    currentNull = Null2;
                    probMassAccounted = probMassAccounted + Alpha_k(k) + Gamma_k(k);
                    if probMassAccounted > 1 + 1e-9
                        warning('Accumulated probability mass > 1 at stage k=%d. Stopping.', k);
                        aThresholds(k+1:end) = 1; gThresholds(k+1:end) = 0; break;
                    end
                end
                stageAlphas = aThresholds; stageGammas = gThresholds;
                obj.thresholdCache(cacheKey) = struct('stageAlphas', stageAlphas, 'stageGammas', stageGammas);
            end
        end
    end
end