classdef ORDCalculator < handle
    % Calculates Magnitude Squared Coherence (MSC) from processed time-domain signals.

    properties (SetAccess = protected)
        fs % Sampling frequency (Hz)
        nfft
        nBins
        latestMSC % Stores MSC for the last calculation: (nBins, K_stages, nChannels)
        latestEpochs
        latestK
        latestM
    end

    methods
        function obj = ORDCalculator(fs)
            arguments
                fs (1,1) {mustBeNumeric, mustBePositive}
            end
            obj.fs = fs;
            obj.nfft = fs;
            obj.nBins = floor(obj.nfft / 2) + 1;
        end

        function [msc, epochs, M, K] = calculateMSC(obj, processedSignals, options)
            arguments
                obj
                processedSignals (:,:,:) {mustBeNumeric}
                options.method {mustBeMember(options.method, {'fixedStep', 'fixedK'})} = 'fixedStep'
                options.startWindow (1,1) {mustBeInteger, mustBePositive} = 1
                options.stepSize (1,1) {mustBeInteger, mustBePositive} = 24
                options.K (1,1) {mustBeInteger, mustBePositive} = 5
                options.maxWindow (1,1) {mustBeInteger, mustBePositive} = size(processedSignals, 2)
            end

            if size(processedSignals, 1) ~= obj.fs
                warning('Input signal Fs (%d samples/window) does not match ORDCalculator Fs (%d).', size(processedSignals, 1), obj.fs);
            end
            if options.startWindow >= options.maxWindow
                warning('startWindow (%d) >= maxWindow (%d). Cannot calculate MSC.', options.startWindow, options.maxWindow);
                msc = []; epochs = []; M = 0; K = 0; obj.latestMSC = []; return;
            end

            [epochs, M, K] = obj.defineEpochs(options.method, options.startWindow, options.maxWindow, options.stepSize, options.K);

            if K < 1 || M < 2
                warning('Epoch definition invalid (K=%d, M=%d). Cannot calculate MSC.', K, M);
                msc = []; obj.latestMSC = []; obj.latestEpochs = epochs; obj.latestK = K; obj.latestM = M; return;
            end

            nChannels = size(processedSignals, 3);
            msc = zeros(obj.nBins, K, nChannels);

            for ch = 1:nChannels
                channelData = squeeze(processedSignals(:, :, ch));
                for k_stage = 1:K
                    stageStartWin = epochs(k_stage, 1);
                    stageEndWin = epochs(k_stage, 2);
                    currentStageSignals = channelData(:, stageStartWin : stageEndWin);

                    Y_stage = fft(currentStageSignals, obj.nfft, 1);
                    Y_stage_pos = Y_stage(1:obj.nBins, :);

                    msc_num = abs(sum(Y_stage_pos, 2)).^2;
                    msc_den = M * sum(abs(Y_stage_pos).^2, 2);

                    valid_den = msc_den > eps;
                    msc(valid_den, k_stage, ch) = msc_num(valid_den) ./ msc_den(valid_den);
                end
            end

            msc(isnan(msc)) = 0;
            msc(isinf(msc)) = 0;

            obj.latestMSC = msc;
            obj.latestEpochs = epochs;
            obj.latestK = K;
            obj.latestM = M;
        end

        function [epochs, M, K] = defineEpochs(obj, method, startWindow, maxWindow, stepSize, K_target)
             arguments
                obj, method, startWindow, maxWindow, stepSize, K_target
            end
            totalAvailableWindows = maxWindow - startWindow + 1;
            epochs = []; M = 0; K = 0;

            switch method
                case 'fixedStep'
                    M = stepSize;
                    if M < 2, warning('stepSize (%d) must be >= 2.', M); return; end
                    K = floor(totalAvailableWindows / M);
                    if K == 0, warning('Not enough windows (%d) for stepSize=%d.', totalAvailableWindows, M); return; end
                    epochs = zeros(K, 2);
                    for k_stage = 1:K
                        epochStart = startWindow + (k_stage-1) * M;
                        epochs(k_stage, :) = [epochStart, epochStart + M - 1];
                    end
                case 'fixedK'
                    K = K_target;
                    if K < 1, warning('Requested K (%d) must be >= 1.', K_target); return; end
                    if K > totalAvailableWindows, warning('Requested K (%d) > available windows (%d). Reducing K.', K, totalAvailableWindows); K = totalAvailableWindows; end

                    boundaries = round(linspace(startWindow-1, maxWindow, K + 1));
                    epochs = zeros(K, 2);
                    actualM = zeros(K, 1);
                    for k_stage = 1:K
                        epochs(k_stage, 1) = boundaries(k_stage) + 1;
                        epochs(k_stage, 2) = boundaries(k_stage + 1);
                        actualM(k_stage) = epochs(k_stage, 2) - epochs(k_stage, 1) + 1;
                    end
                    if any(actualM < 2), warning('Epoch definition for fixedK=%d resulted in stages with M < 2.', K); end
                    M = floor(mean(actualM));
                otherwise, error('Unknown epoch definition method: %s', method);
            end
        end

         function [msc, epochs, M, K] = getLatestResult(obj)
             msc = obj.latestMSC; epochs = obj.latestEpochs; M = obj.latestM; K = obj.latestK;
         end
    end
end