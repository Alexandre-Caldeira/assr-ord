% --- FILE START: CGST_BetaStrategy.m ---
classdef CGST_BetaStrategy < TestStrategy
    % Implements the Chesnaye Group Sequential Test (CGST) strategy using
    % thresholds derived from the Beta(1, M-1) distribution under H0 (MSC null distribution).
    % Requires the Statistics and Machine Learning Toolbox for betapdf.

    properties
        RequiresFutilityThreshold = true % This strategy calculates both efficacy and futility bounds.

        % Cache for thresholds to avoid re-computation. Key: 'K<k>_M<m>_a<a>'.
        thresholdCache = containers.Map('KeyType', 'char', 'ValueType', 'any');

        Name = 'CGST_Beta' % Strategy identifier
    end

    methods
        function obj = CGST_BetaStrategy()
            % Constructor - initializes the threshold cache.
            obj.thresholdCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
            fprintf('CGST_BetaStrategy: Initialized. Requires Statistics Toolbox for betapdf.\n');
        end

        function thresholds = calculateThresholds(obj, strategyParams, analysisInfo)
            % Calculates CGST efficacy (alpha) and futility (gamma) thresholds.
            % Args:
            %   strategyParams: Struct must contain .alpha (desired overall type I error).
            %   analysisInfo: Struct must contain .K (number of stages) and .M (windows per stage).
            %                 .M is assumed constant across stages for standard CGST.
            % Returns:
            %   thresholds: Struct with .efficacy and .futility threshold vectors (length K).

            % --- Input Validation ---
            if ~isfield(strategyParams, 'alpha') || ~isscalar(strategyParams.alpha) || strategyParams.alpha <= 0 || strategyParams.alpha >= 1
                error('CGST_BetaStrategy: strategyParams must include a scalar alpha between 0 and 1.');
            end
            if ~isfield(analysisInfo, 'K') || ~isscalar(analysisInfo.K) || analysisInfo.K < 1 || floor(analysisInfo.K) ~= analysisInfo.K
                 error('CGST_BetaStrategy: analysisInfo must include a positive integer K (number of stages).');
            end
            if ~isfield(analysisInfo, 'M') || ~isscalar(analysisInfo.M) || analysisInfo.M <= 1 || floor(analysisInfo.M) ~= analysisInfo.M
                % Allow vector M but warn and use M(1) - standard CGST assumes constant M
                 if isfield(analysisInfo, 'M') && isvector(analysisInfo.M) && numel(analysisInfo.M) == analysisInfo.K && all(analysisInfo.M > 1)
                     warning('CGST_BetaStrategy: analysisInfo.M is a vector. Standard CGST assumes constant M. Using M=%d from the first stage for threshold calculation.', analysisInfo.M(1));
                     M = analysisInfo.M(1);
                 else
                    error('CGST_BetaStrategy: analysisInfo must include a positive integer M > 1 (windows per stage). Vector M is not standard.');
                 end
            else
                M = analysisInfo.M; % M is scalar and valid
            end
            K = analysisInfo.K;
            alpha = strategyParams.alpha;

            % --- Check Cache ---
            cacheKey = sprintf('K%d_M%d_a%.6f', K, M, alpha); % Use more precision for alpha key
            if isKey(obj.thresholdCache, cacheKey)
                cachedData = obj.thresholdCache(cacheKey);
                thresholds.efficacy = cachedData.a;
                thresholds.futility = cachedData.g;
                % fprintf('CGST_BetaStrategy: Loaded thresholds from cache for %s\n', cacheKey);
                return;
            end
            fprintf('CGST_BetaStrategy: Calculating thresholds for %s...\n', cacheKey);

            % --- Parameters for Calculation ---
            aThresholds = zeros(1, K); % Efficacy thresholds
            gThresholds = zeros(1, K); % Futility thresholds
            Alpha_k = ones(1, K) * (alpha / K); % Equal alpha spending per stage
            Gamma_k = ((1 - alpha) / K) * ones(1, K); % Equal gamma spending per stage

            Resolution = 1e5; % Resolution for PDF approximation (higher is more accurate but slower)
            Xvalues = linspace(0, 1, Resolution + 1); % Values from 0 to 1

            % --- Base Null Distribution (Beta(1, M-1)) ---
            try
                Null = betapdf(Xvalues, 1, M - 1);
                 % Handle potential NaN/Inf at boundaries (esp. if M is close to 1, though M>1 enforced)
                 Null(isnan(Null) | isinf(Null)) = 0;
                 if sum(Null) < eps
                      error('CGST_BetaStrategy: Initial Beta(1, %d) PDF sum is near zero. Check parameters.', M-1);
                 end
                 Null = Null / sum(Null); % Normalize
                 Chi2_Norm = Null; % Base distribution for convolution
            catch ME
                if contains(ME.identifier, 'stats') || contains(ME.message, 'Statistics')
                    error('CGST_BetaStrategy: betapdf function from Statistics Toolbox is required but not found or errored: %s', ME.message);
                else
                    rethrow(ME); % Rethrow other errors
                end
            end

            % --- Stage 1 Thresholds (Closed-form) ---
            aThresholds(1) = 1 - Alpha_k(1)^(1 / (M - 1));
            gThresholds(1) = 1 - (1 - Gamma_k(1))^(1 / (M - 1));
            % Clamp thresholds to [0, 1] for safety
            aThresholds(1) = max(0, min(1, aThresholds(1)));
            gThresholds(1) = max(0, min(1, gThresholds(1)));

            % Convert thresholds to indices for truncation (+1 for 1-based index)
            TruncInd_Ra = floor(aThresholds(1) * Resolution) + 1;
            TruncInd_Rg = floor(gThresholds(1) * Resolution) + 1;
            % Clamp indices to valid range [1, Resolution+1]
            TruncInd_Ra = min(TruncInd_Ra, Resolution + 1);
            TruncInd_Rg = max(TruncInd_Rg, 1);

            % --- Iterative Threshold Calculation (Stages 2 to K) ---
            currentNull = Null; % PDF of the cumulative statistic at the start of the stage
            prob_mass_carried_forward = 1.0;

            for k = 2:K
                % Truncate the distribution based on *previous* stage's boundaries
                NullTrunc = currentNull;
                if TruncInd_Ra <= numel(NullTrunc); NullTrunc(TruncInd_Ra:end) = 0; end
                if TruncInd_Rg >= 1; NullTrunc(1:TruncInd_Rg) = 0; end

                % Check if truncation eliminated all probability mass
                sumTrunc = sum(NullTrunc);
                if sumTrunc < eps
                     warning('CGST_BetaStrategy: Null distribution truncated to zero before stage %d convolution. Thresholds for stage >= %d set to NaN. Check alpha/gamma allocation or M.', k, k);
                     aThresholds(k:end) = NaN;
                     gThresholds(k:end) = NaN;
                     break; % Stop calculation
                end
                 NullTrunc = NullTrunc / sumTrunc; % Renormalize truncated part

                % Convolve with the base distribution
                % Length of result is length(A)+length(B)-1
                Null2 = conv(Chi2_Norm, NullTrunc, 'full');
                if sum(Null2) < eps
                     warning('CGST_BetaStrategy: Convolution resulted in zero sum at stage %d. Thresholds for stage >= %d set to NaN.', k, k);
                     aThresholds(k:end) = NaN;
                     gThresholds(k:end) = NaN;
                     break; % Stop calculation
                end
                Null2 = Null2 / sum(Null2); % Normalize convolution result


                % Adjust total probability mass remaining for this stage's calculation
                prob_mass_entering_stage_k = prob_mass_carried_forward * sumTrunc; % Mass that didn't stop at k-1
                if prob_mass_entering_stage_k < eps
                     warning('CGST_BetaStrategy: Probability mass entering stage %d is near zero. Thresholds for stage >= %d set to NaN.', k, k);
                     aThresholds(k:end) = NaN;
                     gThresholds(k:end) = NaN;
                     break;
                end

                % Find new thresholds based on the PDF of the statistic *at stage k* (Null2)
                cumulativeNull2 = cumsum(Null2);
                % Target percentiles within the distribution of those *continuing* to stage k
                targetAlphaPercentile = 1 - Alpha_k(k) / prob_mass_entering_stage_k;
                targetGammaPercentile = Gamma_k(k) / prob_mass_entering_stage_k;

                 % Clamp target percentiles to [0, 1] to avoid issues if spending exceeds remaining mass
                 targetAlphaPercentile = max(0, min(1, targetAlphaPercentile));
                 targetGammaPercentile = max(0, min(1, targetGammaPercentile));


                % Find indices corresponding to these percentiles
                idx_a = find(cumulativeNull2 >= targetAlphaPercentile, 1, 'first');
                idx_g = find(cumulativeNull2 >= targetGammaPercentile, 1, 'first');

                % Handle cases where find returns empty (e.g., percentile is 0 or 1 exactly)
                convResolution = numel(Null2) - 1; % Effective resolution of the convolved PDF scale [0, convResolution]
                if isempty(idx_a); idx_a = numel(Null2); end % If percentile is 1 (or > max), set to end
                if isempty(idx_g); idx_g = 1; end % If percentile is 0 (or < min), set to start

                % Convert indices back to thresholds [0, 1] scale
                % The convolved scale is roughly [0, 2], so map index to [0, 1]
                % A value at index `idx` corresponds to roughly `(idx-1)/convResolution` on the original scale
                % Need to map the convolved scale back. The effective range of the sum of k vars each ~[0,1] is ~[0,k].
                % Using the helper function findIndex logic from ORDTester might be more robust here if available,
                % or carefully map the convolved index back to the [0,1] MSC scale.
                % Approximation: Assume convolution stretches the scale linearly.
                aThresholds(k) = (idx_a - 1) / convResolution; % This assumes the convolved PDF maps linearly to [0,1] scale, which is an approximation.
                gThresholds(k) = (idx_g - 1) / convResolution;

                 % Clamp thresholds and update indices for next iteration's truncation
                 aThresholds(k) = max(0, min(1, aThresholds(k)));
                 gThresholds(k) = max(0, min(1, gThresholds(k)));
                 TruncInd_Ra = floor(aThresholds(k) * convResolution) + 1;
                 TruncInd_Rg = floor(gThresholds(k) * convResolution) + 1;
                  % Clamp indices
                 TruncInd_Ra = min(TruncInd_Ra, convResolution + 1);
                 TruncInd_Rg = max(TruncInd_Rg, 1);

                 currentNull = Null2; % Update the distribution for the next stage
                 prob_mass_carried_forward = prob_mass_entering_stage_k * (1 - targetGammaPercentile - (1-targetAlphaPercentile)); % Update remaining mass


            end % End loop for stages k=2:K

            thresholds.efficacy = aThresholds;
            thresholds.futility = gThresholds;

            % Cache the result
            obj.thresholdCache(cacheKey) = struct('a', aThresholds, 'g', gThresholds);
            fprintf('CGST_BetaStrategy: Threshold calculation complete and cached for %s.\n', cacheKey);
        end


        function [decisionSequence, stoppingStage] = runTest(obj, statSequence, thresholds)
            % Runs the CGST test using cumulative MSC values.
            % Args:
            %   statSequence: MSC values (bins x stages x channels).
            %   thresholds: Struct with .efficacy and .futility vectors (length K).
            % Returns:
            %   decisionSequence: Decisions (1=Detect, -1=Futility, 0=Continue) (bins x stages x channels).
            %   stoppingStage: Stage number of decision (bins x channels).

            if ~isfield(thresholds, 'efficacy') || ~isfield(thresholds, 'futility')
                error('CGST_BetaStrategy: Threshold struct missing efficacy or futility fields.');
            end
            if isempty(statSequence); decisionSequence=[]; stoppingStage=[]; return; end

            [nBins, K, nChannels] = size(statSequence);
            aThresh = thresholds.efficacy(:)'; % Ensure row vector
            gThresh = thresholds.futility(:)'; % Ensure row vector

            if length(aThresh) ~= K || length(gThresh) ~= K
                error('CGST_BetaStrategy: Threshold vector length (%d) does not match number of stages (%d)', length(aThresh), K);
            end

            decisionSequence = zeros(nBins, K, nChannels); % 0 = continue initially
            stoppingStage = ones(nBins, nChannels) * K;    % Default: stops at last stage K

            % --- Calculate Cumulative Statistic ---
            % Assumes statSequence contains the *per-stage* statistic (e.g., MSC from stage k alone)
            cumulativeStat = cumsum(statSequence, 2); % Accumulate across stages (dim 2)

            % --- Apply Thresholds Stage by Stage ---
            for chan = 1:nChannels
                for freqBin = 1:nBins
                    alreadyStopped = false;
                    for k = 1:K
                        if alreadyStopped
                            % Carry forward the decision from the stopping stage
                            decisionSequence(freqBin, k, chan) = decisionSequence(freqBin, k-1, chan);
                            continue;
                        end

                        currentCumulative = cumulativeStat(freqBin, k, chan);

                        % Check for NaN input -> NaN decision output
                        if isnan(currentCumulative)
                             decisionSequence(freqBin, k:K, chan) = NaN; % Propagate NaN decision
                             stoppingStage(freqBin, chan) = k; % Mark stop stage
                             alreadyStopped = true;
                             continue;
                        end

                        % Check Efficacy Threshold
                        if currentCumulative >= aThresh(k)
                            decisionSequence(freqBin, k:K, chan) = 1; % Detect from this stage onwards
                            stoppingStage(freqBin, chan) = k;
                            alreadyStopped = true;
                            continue; % Move to next freq/chan
                        end

                        % Check Futility Threshold
                        if currentCumulative < gThresh(k)
                            decisionSequence(freqBin, k:K, chan) = -1; % Futility from this stage onwards
                            stoppingStage(freqBin, chan) = k;
                            alreadyStopped = true;
                            continue; % Move to next freq/chan
                        end

                        % If neither threshold crossed, decision remains 0 (Continue)
                        % stoppingStage remains K until a decision is made or loop ends
                    end % stage loop (k)
                end % frequency loop (freqBin)
            end % channel loop (chan)
        end % runTest method
    end % methods
end % classdef