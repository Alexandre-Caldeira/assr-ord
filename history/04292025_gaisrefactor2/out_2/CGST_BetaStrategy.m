% --- FILE START: CGST_BetaStrategy.m ---
classdef CGST_BetaStrategy < TestStrategy
    % Implements the Chesnaye Group Sequential Test (CGST) strategy using
    % thresholds derived from the Beta(1, M-1) distribution under H0 (MSC null distribution).
    % Requires the Statistics and Machine Learning Toolbox for betapdf.

    properties
        % Implementing abstract properties:
        RequiresFutilityThreshold = true
        Name = 'CGST_Beta'

        % Specific properties for this class:
        thresholdCache % Cache for thresholds. Key: 'K<k>_M<m>_a<a>'.
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
                 if isfield(analysisInfo, 'M') && isvector(analysisInfo.M) && ~isempty(analysisInfo.M) && all(analysisInfo.M > 1)
                     if numel(unique(analysisInfo.M)) > 1
                        warning('CGST_BetaStrategy: analysisInfo.M is a vector with varying values. Standard CGST assumes constant M. Using M=%d from the first stage.', analysisInfo.M(1));
                     end
                     M = analysisInfo.M(1); % Use first M value
                 else
                    error('CGST_BetaStrategy: analysisInfo must include a positive integer M > 1 (windows per stage). Check vector M if provided.');
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
                return;
            end
            fprintf('CGST_BetaStrategy: Calculating thresholds for %s...\n', cacheKey);

            % --- Parameters for Calculation ---
            aThresholds = zeros(1, K); gThresholds = zeros(1, K);
            Alpha_k = ones(1, K) * (alpha / K);
            Gamma_k = ((1 - alpha) / K) * ones(1, K);

            Resolution = 1e5;
            Xvalues = linspace(0, 1, Resolution + 1);

            % --- Base Null Distribution (Beta(1, M-1)) ---
            try
                Null = betapdf(Xvalues, 1, M - 1);
                 Null(isnan(Null) | isinf(Null)) = 0;
                 if sum(Null) < eps; error('CGST_BetaStrategy: Initial Beta(1, %d) PDF sum is near zero.', M-1); end
                 Null = Null / sum(Null); % Normalize
                 Chi2_Norm = Null;
            catch ME
                if contains(ME.identifier, 'stats') || contains(ME.message, 'Statistics')
                    error('CGST_BetaStrategy: betapdf function from Statistics Toolbox required: %s', ME.message);
                else; rethrow(ME); end
            end

            % --- Stage 1 Thresholds ---
            aThresholds(1) = max(0, min(1, 1 - Alpha_k(1)^(1 / (M - 1))));
            gThresholds(1) = max(0, min(1, 1 - (1 - Gamma_k(1))^(1 / (M - 1))));

            TruncInd_Ra = floor(aThresholds(1) * Resolution) + 1;
            TruncInd_Rg = floor(gThresholds(1) * Resolution) + 1;
            TruncInd_Ra = min(TruncInd_Ra, Resolution + 1);
            TruncInd_Rg = max(TruncInd_Rg, 1);

            % --- Iterative Threshold Calculation (Stages 2 to K) ---
            currentNull = Null;
            prob_mass_carried_forward = 1.0;

            for k = 2:K
                NullTrunc = currentNull;
                if TruncInd_Ra <= numel(NullTrunc); NullTrunc(TruncInd_Ra:end) = 0; end
                if TruncInd_Rg >= 1; NullTrunc(1:TruncInd_Rg) = 0; end

                sumTrunc = sum(NullTrunc);
                if sumTrunc < eps
                     warning('CGST_BetaStrategy: Null distribution truncated to zero before stage %d. Thresholds >= %d set to NaN.', k, k);
                     aThresholds(k:end) = NaN; gThresholds(k:end) = NaN; break;
                end
                 NullTrunc = NullTrunc / sumTrunc; % Renormalize

                Null2 = conv(Chi2_Norm, NullTrunc, 'full');
                if sum(Null2) < eps
                     warning('CGST_BetaStrategy: Convolution resulted in zero sum at stage %d. Thresholds >= %d set to NaN.', k, k);
                     aThresholds(k:end) = NaN; gThresholds(k:end) = NaN; break;
                end
                Null2 = Null2 / sum(Null2);

                prob_mass_entering_stage_k = prob_mass_carried_forward * sumTrunc;
                if prob_mass_entering_stage_k < eps
                     warning('CGST_BetaStrategy: Probability mass entering stage %d near zero. Thresholds >= %d set to NaN.', k, k);
                     aThresholds(k:end) = NaN; gThresholds(k:end) = NaN; break;
                end

                cumulativeNull2 = cumsum(Null2);
                targetAlphaPercentile = max(0, min(1, 1 - Alpha_k(k) / prob_mass_entering_stage_k));
                targetGammaPercentile = max(0, min(1, Gamma_k(k) / prob_mass_entering_stage_k));

                idx_a = find(cumulativeNull2 >= targetAlphaPercentile, 1, 'first');
                idx_g = find(cumulativeNull2 >= targetGammaPercentile, 1, 'first');

                convResolution = numel(Null2) - 1;
                if isempty(idx_a); idx_a = numel(Null2); end
                if isempty(idx_g); idx_g = 1; end

                % Use ORDTester.findIndex logic if available for potentially more accuracy
                % Otherwise, use linear mapping approximation
                if exist('ORDTester', 'class') && ismethod('ORDTester', 'findIndex')
                    try
                        idx_a_alt = ORDTester.findIndex(Null2, sum(Null2) - Alpha_k(k)/prob_mass_entering_stage_k); % Find index for efficacy
                        idx_g_alt = ORDTester.findIndex(Null2, Gamma_k(k)/prob_mass_entering_stage_k, 1); % Find index for futility (GetFutil=1)
                        aThresholds(k) = (idx_a_alt - 1) / convResolution; % Adjust index and scale
                        gThresholds(k) = (idx_g_alt - 1) / convResolution;
                    catch ME_find
                         warning('CGST_BetaStrategy: Failed to use ORDTester.findIndex, using linear approx. Error: %s', ME_find.message);
                         aThresholds(k) = (idx_a - 1) / convResolution;
                         gThresholds(k) = (idx_g - 1) / convResolution;
                    end
                else
                    aThresholds(k) = (idx_a - 1) / convResolution;
                    gThresholds(k) = (idx_g - 1) / convResolution;
                end


                 aThresholds(k) = max(0, min(1, aThresholds(k)));
                 gThresholds(k) = max(0, min(1, gThresholds(k)));
                 TruncInd_Ra = floor(aThresholds(k) * convResolution) + 1;
                 TruncInd_Rg = floor(gThresholds(k) * convResolution) + 1;
                 TruncInd_Ra = min(TruncInd_Ra, convResolution + 1);
                 TruncInd_Rg = max(TruncInd_Rg, 1);

                 currentNull = Null2;
                 prob_mass_carried_forward = prob_mass_entering_stage_k * (targetAlphaPercentile - targetGammaPercentile); % Remaining probability

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
                error('CGST_BetaStrategy: Threshold vector length (%d, %d) does not match number of stages (%d)', length(aThresh), length(gThresh), K);
            end

            decisionSequence = zeros(nBins, K, nChannels); % 0 = continue initially
            stoppingStage = ones(nBins, nChannels) * K;    % Default: stops at last stage K

            cumulativeStat = cumsum(statSequence, 2); % Accumulate across stages (dim 2)

            for chan = 1:nChannels
                for freqBin = 1:nBins
                    alreadyStopped = false;
                    for k = 1:K
                        if alreadyStopped
                            decisionSequence(freqBin, k, chan) = decisionSequence(freqBin, k-1, chan);
                            continue;
                        end

                        currentCumulative = cumulativeStat(freqBin, k, chan);

                        if isnan(currentCumulative) || isnan(aThresh(k)) || isnan(gThresh(k))
                             decisionSequence(freqBin, k:K, chan) = NaN; % Propagate NaN decision if input or threshold is NaN
                             stoppingStage(freqBin, chan) = k;
                             alreadyStopped = true;
                             continue;
                        end

                        % Check Efficacy Threshold
                        if currentCumulative >= aThresh(k)
                            decisionSequence(freqBin, k:K, chan) = 1;
                            stoppingStage(freqBin, chan) = k;
                            alreadyStopped = true;
                            continue;
                        end

                        % Check Futility Threshold
                        if currentCumulative < gThresh(k)
                            decisionSequence(freqBin, k:K, chan) = -1;
                            stoppingStage(freqBin, chan) = k;
                            alreadyStopped = true;
                            continue;
                        end
                    end % stage loop (k)
                end % frequency loop (freqBin)
            end % channel loop (chan)
        end % runTest method
    end % methods
end % classdef