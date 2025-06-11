classdef ORDTester < handle
    % Applies statistical tests (e.g., Beta CGST) to MSC data for objective response detection.
    % Calculates thresholds and determines detection outcomes (TP, TN, FP, FN).

    properties (SetAccess = protected)
        desiredAlpha = 0.05; % Overall Type I error rate (alpha) for the test sequence

        % Threshold Cache (stores calculated thresholds to avoid recomputation)
        % Key: M (windows per stage), K (number of stages)
        % Value: Struct with fields 'stageAlphas' and 'stageGammas'
        thresholdCache % Use a Map for potentially better lookup than cell array

        % Frequency Information (needed to classify TP/FN vs FP/TN)
        signalFrequencies
        noiseFrequencies
        allTestFrequencies
        nSignalFrequencies
        nNoiseFrequencies
        nAllTestFrequencies
        noiseFlagIndex % Index in allTestFrequencies where noise starts
    end

    methods
        function obj = ORDTester(signalFrequencies, noiseFrequencies, options)
            % Constructor for ORDTester
            % signalFrequencies: Vector of signal frequencies
            % noiseFrequencies: Vector of noise frequencies
            % options: Name-Value pairs:
            %   'desiredAlpha': Overall Type I error rate

            arguments % Define fixed args and options block
                signalFrequencies (1,:) {mustBeNumeric , mustBePositive}
                noiseFrequencies (1,:) {mustBeNumeric , mustBePositive} % Allow empty

                % Optional Name-Value arguments
                options.desiredAlpha (1,1) {mustBeNumeric, ....
                            mustBeGreaterThan(options.desiredAlpha,0), ...
                            mustBeLessThan(options.desiredAlpha,1)} = 0.05
            end

            % Assign fixed arguments
            obj.signalFrequencies = sort(signalFrequencies);
            obj.noiseFrequencies = sort(noiseFrequencies);

            % Assign optional arguments
            obj.desiredAlpha = options.desiredAlpha;
            obj.allTestFrequencies = sort([obj.signalFrequencies, obj.noiseFrequencies]);
            
            commonFreqs = intersect(obj.signalFrequencies, obj.noiseFrequencies);
            if ~isempty(commonFreqs)
                error('Signal frequencies and noise frequencies must be distinct. Common: %s', mat2str(commonFreqs));
            end

            obj.nSignalFrequencies = numel(obj.signalFrequencies);
            obj.nNoiseFrequencies = numel(obj.noiseFrequencies);
            obj.nAllTestFrequencies = numel(obj.allTestFrequencies);

            % Find index where noise frequencies begin
            if obj.nNoiseFrequencies > 0
                [~, obj.noiseFlagIndex] = ismember(min(obj.noiseFrequencies), obj.allTestFrequencies);
                 if obj.noiseFlagIndex == 0 % Should not happen if logic is correct
                     error('Could not determine noise frequency start index.');
                 end
            else
                obj.noiseFlagIndex = obj.nAllTestFrequencies + 1; % No noise frequencies
            end


            % Initialize threshold cache
            obj.thresholdCache = containers.Map('KeyType', 'char', 'ValueType', 'any');

            fprintf('[%s] ORDTester initialized. Desired Alpha=%.3f. Signal Freqs: %d, Noise Freqs: %d.\n', ...
                    datetime, obj.desiredAlpha, obj.nSignalFrequencies, obj.nNoiseFrequencies);
        end

        function [decisionResults, stageAlphas, stageGammas] = runBetaCGST(obj, mscData, M, K)
            % Runs the Beta distribution-based CGST on the provided MSC data.
            % mscData: MSC values (nBins, K_stages, nChannels) from ORDCalculator
            % M: Number of windows per stage (used for threshold calculation)
            % K: Number of stages (used for threshold calculation)
            % Returns:
            %   decisionResults: Struct containing TP, TN, FP, FN matrices (nFreqs, K_stages, nChannels)
            %                    indicating outcome *at each stage* for each tested frequency.
            %   stageAlphas: Alpha thresholds used for each stage (1xK vector)
            %   stageGammas: Gamma thresholds used for each stage (1xK vector)

            arguments
                obj
                mscData (:,:,:) {mustBeNumeric} % nBins x K x nChannels
                M (1,1) {mustBeInteger, mustBePositive}
                K (1,1) {mustBeInteger, mustBePositive}
            end

            if K ~= size(mscData, 2)
                error('Number of stages K (%d) does not match MSC data dimension 2 (%d).', K, size(mscData, 2));
            end
             if M < 2
                error('M (windows per stage) must be >= 2 for Beta CGST thresholds. Got M=%d.', M);
             end

            fprintf('[%s] Running Beta CGST: M=%d, K=%d, Alpha=%.3f...\n', datetime, M, K, obj.desiredAlpha);

            % --- 1. Get or Calculate Thresholds ---
            [stageAlphas, stageGammas] = obj.getBetaCGSTThresholds(M, K);

            % --- 2. Extract MSC at Test Frequencies ---
            % Find indices of test frequencies within the MSC bins (assuming Fs and nfft match)
            % Bins correspond to frequencies 0, df, 2*df, ..., Fs/2
            % df = Fs / nfft. If nfft=Fs, df=1 Hz. Bin index = Freq + 1.
            freqIndices = obj.allTestFrequencies + 1; % +1 for MATLAB 1-based indexing
            maxBinIndex = size(mscData, 1);
            if any(freqIndices > maxBinIndex)
                warning('Some test frequencies (%s) exceed max bin index (%d). Clipping.', ...
                        mat2str(obj.allTestFrequencies(freqIndices > maxBinIndex)), maxBinIndex);
                freqIndices = min(freqIndices, maxBinIndex);
            end

            mscAtTestFreqs = mscData(freqIndices, :, :); % nTestFreqs x K x nChannels

            % --- 3. Apply Stage-Wise Decisions ---
            nChannels = size(mscData, 3);
            TP = zeros(obj.nAllTestFrequencies, K, nChannels);
            TN = zeros(obj.nAllTestFrequencies, K, nChannels);
            FP = zeros(obj.nAllTestFrequencies, K, nChannels);
            FN = zeros(obj.nAllTestFrequencies, K, nChannels);

            % Pre-calculate cumulative sum of MSC across stages for efficiency
            % cumsum operates along dimension 2 (stages)
            cumulativeMSC = cumsum(mscAtTestFreqs, 2);

            for ch = 1:nChannels
                for k_stage = 1:K
                    currentCumulativeMSC = cumulativeMSC(:, k_stage, ch); % Vector for all freqs at this stage

                    % Conditions for decision at stage k
                    detected = currentCumulativeMSC > stageAlphas(k_stage);
                    stopped = currentCumulativeMSC <= stageGammas(k_stage); % Futility/stop

                    % Identify signal and noise frequencies among the test frequencies
                    isSignal = false(obj.nAllTestFrequencies, 1);
                    isSignal(1 : obj.noiseFlagIndex - 1) = true; % Frequencies before noise flag are signal
                    isNoise = ~isSignal;

                    % --- Assign outcomes ---
                    % TP: Signal detected
                    TP(isSignal & detected, k_stage, ch) = 1;
                    % FN: Signal stopped (futility)
                    FN(isSignal & stopped, k_stage, ch) = 1;
                    % FP: Noise detected
                    FP(isNoise & detected, k_stage, ch) = 1;
                    % TN: Noise stopped (futility)
                    TN(isNoise & stopped, k_stage, ch) = 1;

                    % --- Handle ongoing tests ---
                    % If a frequency was detected or stopped at an earlier stage,
                    % it shouldn't contribute to outcomes in later stages.
                    % We need to track the *first* stage a decision is made.

                    % Let's redefine: TP/TN/FP/FN(f, k) = 1 if decision for freq f is made *exactly* at stage k.
                    % Reset current stage decisions
                    TP(:, k_stage, ch) = 0; FN(:, k_stage, ch) = 0;
                    FP(:, k_stage, ch) = 0; TN(:, k_stage, ch) = 0;

                    % Find frequencies decided *before* this stage
                    if k_stage > 1
                        decided_earlier = any(TP(:, 1:k_stage-1, ch), 2) | ...
                                          any(FN(:, 1:k_stage-1, ch), 2) | ...
                                          any(FP(:, 1:k_stage-1, ch), 2) | ...
                                          any(TN(:, 1:k_stage-1, ch), 2);
                    else
                        decided_earlier = false(obj.nAllTestFrequencies, 1);
                    end

                    % Consider only frequencies not decided earlier
                    undecided = ~decided_earlier;

                    % Assign outcomes for undecided frequencies at this stage k
                    TP(undecided & isSignal & detected, k_stage, ch) = 1;
                    FN(undecided & isSignal & stopped, k_stage, ch) = 1;
                    FP(undecided & isNoise & detected, k_stage, ch) = 1;
                    TN(undecided & isNoise & stopped, k_stage, ch) = 1;

                     % Check for frequencies neither detected nor stopped by final stage
                     if k_stage == K
                         final_undecided = undecided & ~detected & ~stopped;
                         if any(final_undecided)
                             % How to classify these? Often considered non-detections.
                             % Assign based on signal/noise status.
                             FN(final_undecided & isSignal, k_stage, ch) = 1; % Signal present but not detected -> FN
                             TN(final_undecided & isNoise, k_stage, ch) = 1;  % Noise present and not detected -> TN
                             % fprintf('  Warning: %d frequencies undecided at final stage K=%d for channel %d. Classified as FN/TN.\n', ...
                             %        sum(final_undecided), K, ch);
                         end
                     end
                end % end stage loop
            end % end channel loop

            decisionResults = struct('TP', TP, 'TN', TN, 'FP', FP, 'FN', FN, ...
                                     'allTestFrequencies', obj.allTestFrequencies, ...
                                     'signalFrequencies', obj.signalFrequencies, ...
                                     'noiseFrequencies', obj.noiseFrequencies);
            fprintf('[%s] Beta CGST decisions calculated.\n', datetime);
        end


        function [stageAlphas, stageGammas] = getBetaCGSTThresholds(obj, M, K)
            % Calculates or retrieves cached Beta CGST thresholds for given M and K.
            arguments
                obj
                M (1,1) {mustBeInteger, mustBeGreaterThanOrEqual(M,2)} % Windows per stage
                K (1,1) {mustBeInteger, mustBePositive} % Number of stages
            end

            cacheKey = sprintf('M%dK%d', M, K);

            if obj.thresholdCache.isKey(cacheKey)
                % Retrieve from cache
                cachedData = obj.thresholdCache(cacheKey);
                stageAlphas = cachedData.stageAlphas;
                stageGammas = cachedData.stageGammas;
                % fprintf(' (Cache Hit) '); % Optional: for debugging cache usage
            else
                % Calculate thresholds
                fprintf('[%s] Calculating Beta CGST thresholds for M=%d, K=%d...', datetime, M, K);
                alpha = obj.desiredAlpha;
                Alpha_k = ones(1, K) * (alpha / K); % Equal alpha allocation per stage

                % Gamma allocation (futility): Equal allocation of remaining probability
                % This is the 'standard' approach, matching the original code's
                % 'single_exam_beta_cgst_threshold' logic.
                Gamma_k = ((1 - alpha) / K) .* ones(1, K);

                % --- Calculation based on original 'single_exam_beta_cgst_threshold' ---
                % High resolution for numerical integration/convolution approximation
                Resolution = 1e5; % (1/0.00001)
                Xvalues = 0 : (1/Resolution) : 1;

                % Null distribution for a single stage (Beta(1, M-1) for MSC under H0)
                % Handle M=2 case where Beta(1,1) is uniform
                 if M == 2
                    Null = ones(1, length(Xvalues)); % Uniform distribution
                 else
                    Null = betapdf(Xvalues, 1, M - 1);
                 end
                 Null = Null / sum(Null); % Normalize to PDF

                Chi2_Norm = Null; % Distribution of a single stage's statistic under H0

                aThresholds = zeros(1, K);
                gThresholds = zeros(1, K);

                % Stage k=1 thresholds
                aThresholds(1) = 1 - Alpha_k(1).^(1 ./ (M - 1)); % Exact for Beta(1, M-1) CDF inverse
                gThresholds(1) = 1 - (1 - Gamma_k(1)).^(1 ./ (M - 1)); % Exact for Beta(1, M-1) CDF inverse

                % --- Iterative calculation for k = 2 to K using convolution ---
                currentNull = Null; % PDF of sum statistic up to stage k-1
                probMassAccounted = Alpha_k(1) + Gamma_k(1);

                for k = 2:K
                    % Truncate the previous null distribution based on stage k-1 thresholds
                    TruncInd_Ra_prev = round(aThresholds(k-1) * Resolution) + 1; % +1 for index
                    TruncInd_Rg_prev = round(gThresholds(k-1) * Resolution) + 1;

                    NullTrunc = currentNull;
                    if TruncInd_Ra_prev <= length(NullTrunc)
                        NullTrunc(TruncInd_Ra_prev:end) = 0; % Zero out mass above alpha threshold
                    end
                     if TruncInd_Rg_prev >= 1
                        NullTrunc(1:TruncInd_Rg_prev) = 0;   % Zero out mass below gamma threshold
                     end
                     
                     % Renormalize truncated distribution (optional but good practice)
                     sumTrunc = sum(NullTrunc);
                     if sumTrunc > 1e-9 % Avoid division by zero
                        NullTrunc = NullTrunc / sumTrunc;
                     else
                         warning('Null distribution truncated to near zero at stage k=%d. Thresholds might be unreliable.', k);
                         % Assign extreme thresholds if distribution collapses
                         aThresholds(k:end) = 1; 
                         gThresholds(k:end) = 0;
                         break; % Stop calculation
                     end


                    % Convolve with the single-stage null distribution
                    % This approximates the distribution of the sum of MSC up to stage k
                    convLength = length(currentNull) + length(Chi2_Norm) - 1;
                    Null2_unscaled = conv(NullTrunc, Chi2_Norm);
                    
                    % Ensure length matches expected resolution for indexing
                    if length(Null2_unscaled) > Resolution+1
                        Null2_unscaled = Null2_unscaled(1:Resolution+1);
                    elseif length(Null2_unscaled) < Resolution+1
                        Null2_unscaled(end+1:Resolution+1) = 0;
                    end

                    % Normalize the convoluted distribution to account for the probability mass
                    % already removed by decisions in previous stages (alpha_1..k-1, gamma_1..k-1)
                    remainingProbMass = 1 - probMassAccounted;
                    currentSum = sum(Null2_unscaled);
                     if currentSum > 1e-9 && remainingProbMass > 1e-9
                        Null2 = Null2_unscaled * (remainingProbMass / currentSum);
                     else
                         warning('Normalization factor near zero at stage k=%d. Thresholds might be unreliable.', k);
                         aThresholds(k:end) = 1; 
                         gThresholds(k:end) = 0;
                         break; % Stop calculation
                     end


                    % Find thresholds for stage k based on the convoluted distribution (Null2)
                    % Alpha threshold: Find index where cumulative sum from *right* equals Alpha_k(k)
                    targetAlphaMass = Alpha_k(k);
                    cumulativeSumRight = cumsum(Null2(end:-1:1)); % Cumulative sum from right
                    idxAlpha = find(cumulativeSumRight >= targetAlphaMass, 1, 'first');
                    if isempty(idxAlpha)
                        aThresholds(k) = 1; % Assign max threshold if mass not reached
                    else
                        % Index needs to be mapped back to original Xvalues index
                        aThresholds(k) = (length(Null2) - idxAlpha) / Resolution;
                    end


                    % Gamma threshold: Find index where cumulative sum from *left* equals Gamma_k(k)
                    targetGammaMass = Gamma_k(k);
                    cumulativeSumLeft = cumsum(Null2);
                    idxGamma = find(cumulativeSumLeft >= targetGammaMass, 1, 'first');
                     if isempty(idxGamma)
                         gThresholds(k) = 0; % Assign min threshold if mass not reached
                     else
                         gThresholds(k) = (idxGamma - 1) / Resolution; % -1 for 0-based value
                     end
                     
                     % Clamp thresholds to [0, 1]
                     aThresholds(k) = min(max(aThresholds(k), 0), 1);
                     gThresholds(k) = min(max(gThresholds(k), 0), 1);

                     % Ensure alpha threshold > gamma threshold
                     if aThresholds(k) <= gThresholds(k)
                         % This can happen due to numerical precision or if alpha/gamma allocations
                         % are very large relative to remaining probability mass.
                         % Adjust slightly or issue warning.
                         % warning('Alpha threshold <= Gamma threshold at stage k=%d. Adjusting.', k);
                         % Simple fix: set gamma slightly lower if they touch/cross
                         if abs(aThresholds(k) - gThresholds(k)) < 1/Resolution
                             gThresholds(k) = max(0, aThresholds(k) - 1/Resolution);
                         end
                         % If still crossed, indicates a problem
                         if aThresholds(k) < gThresholds(k)
                              warning('Threshold calculation error: Alpha < Gamma at stage k=%d.', k);
                              % Could try recalculating with higher resolution or different gamma strategy
                              % For now, set gamma to alpha to avoid impossible region
                              gThresholds(k) = aThresholds(k); 
                         end
                     end


                    % Update for next iteration
                    currentNull = Null2;
                    probMassAccounted = probMassAccounted + Alpha_k(k) + Gamma_k(k);

                end % end k loop

                stageAlphas = aThresholds;
                stageGammas = gThresholds;

                % Store in cache
                obj.thresholdCache(cacheKey) = struct('stageAlphas', stageAlphas, 'stageGammas', stageGammas);
                fprintf(' Done.\n');
            end % end if/else cache
        end % end getBetaCGSTThresholds

    end % end methods

    methods (Static)
        % findIndex function from original code - kept for reference if needed,
        % but threshold finding is now integrated into getBetaCGSTThresholds.
        % If specific behavior of findIndex is required, it can be used there.
        %{
        function Ind = findIndex(PDF1, Goal, GetFutil)
            % Original static method - potentially useful if specific logic needed
            % Note: Assumes PDF1 is normalized, Goal is probability mass.
            arguments
                PDF1 (1,:) {mustBeNumeric}
                Goal (1,1) {mustBeNumeric}
                GetFutil (1,1) {mustBeNumericOrLogical} = 0 % 1 for futility (left tail), 0 for efficacy (right tail)
            end

            L = length(PDF1);
            Resolution = L - 1; % Assuming PDF1 corresponds to 0:1/Res:1

            if GetFutil % Find index from left for Gamma (futility)
                cumulativeSumLeft = cumsum(PDF1);
                Ind = find(cumulativeSumLeft >= Goal, 1, 'first');
                if isempty(Ind)
                    Ind = 0; % Goal not reached, threshold is effectively 0
                end
            else % Find index from right for Alpha (efficacy)
                cumulativeSumRight = cumsum(PDF1(end:-1:1));
                idxFromRight = find(cumulativeSumRight >= Goal, 1, 'first');
                if isempty(idxFromRight)
                    Ind = L; % Goal not reached, threshold is effectively 1 (index L)
                else
                    Ind = L - idxFromRight +1; % Convert index from right to index from left
                end
            end
             % Original findIndex had a binary search logic - the above uses cumsum which is often simpler in MATLAB
             % If the exact binary search behavior is needed, the original static method code can be adapted.
        end
        %}
    end

end % end classdef