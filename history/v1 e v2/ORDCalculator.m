classdef ORDCalculator < handle
    % Calculates Magnitude Squared Coherence (MSC) from processed time-domain signals.
    % Manages epoch definitions for MSC calculation.

    properties (SetAccess = protected)
        fs % Sampling frequency (Hz) - should match data
        nfft % FFT points per window (typically == fs for 1s windows)
        nBins % Number of positive frequency bins (fs/2 + 1)

        latestMSC % Stores MSC for the last calculation: (nBins, K_stages, nChannels)
        latestEpochs % Stores epoch indices used for the last MSC calculation
        latestK % Stores number of stages (K) for the last calculation
        latestM % Stores number of windows per stage (M) for the last calculation
    end

    methods
        function obj = ORDCalculator(fs)
            % Constructor for ORDCalculator
            arguments
                fs (1,1) {mustBeNumeric, mustBePositive} % Sampling frequency of the data
            end
            obj.fs = fs;
            obj.nfft = fs; % Assume 1-second windows
            obj.nBins = floor(obj.nfft / 2) + 1;
            fprintf('[%s] ORDCalculator initialized with Fs = %d Hz (nBins = %d).\n', datetime, obj.fs, obj.nBins);
        end

        function [msc, epochs, M, K] = calculateMSC(obj, processedSignals, p)
            % Calculates MSC for a single exam configuration based on specified epoch strategy.
            % processedSignals: Time-domain data (samples, windows, channels) from PreProcessor.
            % p: Structure with parameters for epoch definition:
            %   p.method: 'fixedStep' (like Zanoteli) or 'fixedK' (like Chesnaye)
            %   p.startWindow: Index of the first window to use (e.g., 1)
            %   p.stepSize: (for 'fixedStep') Number of windows per MSC stage
            %   p.K: (for 'fixedK') Desired number of MSC stages (tests)
            %   p.maxWindow: Index of the last window available in processedSignals

            arguments
                obj
                processedSignals (:,:,:) {mustBeNumeric}
                p.method {mustBeMember(p.method, {'fixedStep', 'fixedK'})} = 'fixedStep'
                p.startWindow (1,1) {mustBeInteger, mustBePositive} = 1
                p.stepSize (1,1) {mustBeInteger, mustBePositive} = 24 % Default for fixedStep
                p.K (1,1) {mustBeInteger, mustBePositive} = 5 % Default for fixedK
                p.maxWindow (1,1) {mustBeInteger, mustBePositive} = size(processedSignals, 2)
            end

            if size(processedSignals, 1) ~= obj.fs
                warning('Input signal Fs (%d samples/window) does not match ORDCalculator Fs (%d). Results may be incorrect.', ...
                        size(processedSignals, 1), obj.fs);
            end

            nAvailableWindows = p.maxWindow - p.startWindow + 1;
            if nAvailableWindows <= 0
                error('startWindow (%d) cannot be greater than maxWindow (%d).', p.startWindow, p.maxWindow);
            end
z
            % --- 1. Determine Epoch Boundaries ---
            [epochs, M, K] = obj.defineEpochs(p.method, p.startWindow, p.maxWindow, p.stepSize, p.K);

            if K < 1
                warning('Epoch definition resulted in K=%d stages. Cannot calculate MSC.', K);
                msc = [];
                obj.latestMSC = [];
                obj.latestEpochs = [];
                obj.latestK = K;
                obj.latestM = M;
                return;
            end
             if M < 2
                warning('Epoch definition resulted in M=%d windows per stage. Need at least 2 for MSC. Cannot calculate MSC.', M);
                msc = [];
                obj.latestMSC = [];
                obj.latestEpochs = epochs;
                obj.latestK = K;
                obj.latestM = M;
                return;
            end

            fprintf('[%s] Calculating MSC: Method=%s, K=%d stages, M=%d windows/stage. Epochs: [%s ... %s].\n', ...
                    datetime, p.method, K, M, num2str(epochs(1,:)), num2str(epochs(end,:)));

            % --- 2. Calculate FFT and MSC ---
            nChannels = size(processedSignals, 3);
            msc = zeros(obj.nBins, K, nChannels);

            for ch = 1:nChannels
                channelData = squeeze(processedSignals(:, :, ch)); % Samples x Windows
                for k_stage = 1:K
                    % Get window indices for this stage
                    stageStartWin = epochs(k_stage, 1);
                    stageEndWin = epochs(k_stage, 2);

                    % Extract time-domain data for this stage's windows
                    % Adjust indices relative to the start of processedSignals
                    currentStageSignals = channelData(:, stageStartWin : stageEndWin); % Samples x M

                    % Calculate FFT for these M windows
                    % fft operates on columns by default
                    Y_stage = fft(currentStageSignals, obj.nfft, 1);
                    Y_stage_pos = Y_stage(1:obj.nBins, :); % Keep only positive frequencies (Bins x M)

                    % Calculate MSC for this stage (Zanotelli formula)
                    msc(:, k_stage, ch) = abs(sum(Y_stage_pos, 2)).^2 ./ (M * sum(abs(Y_stage_pos).^2, 2));
                end
            end
            
            % Handle potential NaNs or Infs resulting from zero power at some frequencies
            msc(isnan(msc)) = 0; 
            msc(isinf(msc)) = 0; % Or handle appropriately, e.g., set to a large value or check denominator

            obj.latestMSC = msc;
            obj.latestEpochs = epochs;
            obj.latestK = K;
            obj.latestM = M;
            fprintf('[%s] MSC calculation complete.\n', datetime);
        end

        function [epochs, M, K] = defineEpochs(obj, method, startWindow, maxWindow, stepSize, K_target)
            % Helper function to define epoch start/end indices based on method.
            % Returns:
            %   epochs: A Kx2 matrix where each row is [start_window_index, end_window_index] for a stage.
            %   M: Number of windows per stage.
            %   K: Actual number of stages created.

            arguments
                obj
                method
                startWindow
                maxWindow
                stepSize
                K_target
            end

            totalAvailableWindows = maxWindow - startWindow + 1;
            epochs = [];
            M = 0;
            K = 0;

            switch method
                case 'fixedStep'
                    % Zanoteli style: Fixed number of windows (M = stepSize) per stage.
                    % K is determined by how many steps fit.
                    M = stepSize;
                    if M > totalAvailableWindows
                        warning('stepSize (%d) > available windows (%d). Cannot create stages.', M, totalAvailableWindows);
                        return;
                    end
                    if M < 2
                         warning('stepSize (%d) must be >= 2 for MSC calculation.', M);
                         return;
                    end

                    K = floor(totalAvailableWindows / M);
                    if K == 0
                       warning('Not enough windows (%d) to form even one stage with stepSize=%d.', totalAvailableWindows, M);
                       return;
                    end

                    epochs = zeros(K, 2);
                    for k_stage = 1:K
                        epochStart = startWindow + (k_stage-1) * M;
                        epochEnd = epochStart + M - 1;
                        epochs(k_stage, :) = [epochStart, epochEnd];
                    end

                case 'fixedK'
                    % Chesnaye style: Fixed number of stages (K = K_target).
                    % M is determined by dividing available windows by K. Windows might be unevenly distributed if needed.
                    K = K_target;
                    if K > totalAvailableWindows
                        warning('Requested K (%d) > available windows (%d). Reducing K.', K, totalAvailableWindows);
                        K = totalAvailableWindows; % Each stage gets 1 window - MSC needs >= 2!
                        % Let calculateMSC handle M<2 warning
                    end
                     if K < 1
                        warning('Requested K (%d) must be >= 1.', K_target);
                        return;
                     end

                    % Use linspace to get approximate boundaries, then round/adjust
                    % This aims for roughly equal M, but ensures integer indices.
                    boundaries = round(linspace(startWindow-1, maxWindow, K + 1)); % K+1 points define K intervals
                    
                    epochs = zeros(K, 2);
                    actualM = zeros(K, 1);
                    for k_stage = 1:K
                        epochs(k_stage, 1) = boundaries(k_stage) + 1;
                        epochs(k_stage, 2) = boundaries(k_stage + 1);
                        actualM(k_stage) = epochs(k_stage, 2) - epochs(k_stage, 1) + 1;
                    end
                    
                    % Check if M is consistent and >= 2
                    if any(actualM < 2)
                         warning('Epoch definition for fixedK=%d resulted in stages with < 2 windows. MSC requires M>=2.', K);
                         % Keep epochs as calculated, let calculateMSC handle it further if needed
                    end
                    M = floor(mean(actualM)); % Report average M, though it might vary slightly
                    if std(actualM) > 1e-6 % Check if M varies
                         fprintf('  Note: Window distribution for fixedK resulted in varying M per stage (approx M=%d).\n', M);
                    end


                otherwise
                    error('Unknown epoch definition method: %s', method);
            end
        end

         function [msc, epochs, M, K] = getLatestResult(obj)
             % Returns the results from the most recent MSC calculation.
             msc = obj.latestMSC;
             epochs = obj.latestEpochs;
             M = obj.latestM;
             K = obj.latestK;
         end
    end
end