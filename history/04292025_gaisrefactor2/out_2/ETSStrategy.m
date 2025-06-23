% --- FILE START: ETSStrategy.m ---
classdef ETSStrategy < TestStrategy
    % Implements the Early Termination Strategy (ETS) based on consecutive
    % detections (NDC) and optionally a futility criterion (Pnv).

    properties
        % Implementing the abstract properties from TestStrategy:
        RequiresFutilityThreshold logical % Set in constructor
        Name char                % Set in constructor

        % Specific properties for this concrete class:
        UseFutility (1,1) logical % Flag to control Pnv logic.
    end

     methods
        function obj = ETSStrategy(useFutility)
            % Constructor: Sets whether to use the futility criterion.
            arguments
                 useFutility logical = false % Default to original ETS (NDC only).
            end
            obj.UseFutility = useFutility;
            obj.RequiresFutilityThreshold = useFutility; % Implement abstract property

            if useFutility
                 obj.Name = 'ETS_Pnv';
            else
                 obj.Name = 'ETS_NDC';
            end
             fprintf('%s: Initialized (UseFutility=%d).\n', obj.Name, obj.UseFutility);
        end % Constructor

        function thresholds = calculateThresholds(obj, strategyParams, analysisInfo)
             % Calculates/retrieves thresholds for the ETS strategy.
             % Args:
             %   strategyParams: Struct must contain .alpha (per-test significance), .NDC (consecutive detections needed).
             %                   If obj.RequiresFutilityThreshold=true, must also contain .VC_not (futility thresholds indexed by M).
             %   analysisInfo: Struct must contain .MM (vector of window counts, M, for each test point/stage).
             % Returns:
             %   thresholds: Struct with .detection (vector), .futility (vector or []), .NDC (scalar).

             % --- Input Validation ---
             if ~isfield(strategyParams, 'alpha') || ~isscalar(strategyParams.alpha) || strategyParams.alpha <= 0 || strategyParams.alpha >= 1
                 error('%s: strategyParams must include a scalar alpha between 0 and 1.', obj.Name);
             end
             if ~isfield(strategyParams, 'NDC') || ~isscalar(strategyParams.NDC) || strategyParams.NDC < 1 || floor(strategyParams.NDC) ~= strategyParams.NDC
                 error('%s: strategyParams must include a positive integer NDC >= 1.', obj.Name);
             end
             if ~isfield(analysisInfo, 'MM') || ~isvector(analysisInfo.MM) || isempty(analysisInfo.MM) || any(analysisInfo.MM < 1)
                  error('%s: analysisInfo must include a non-empty vector MM of positive window counts (M values).', obj.Name);
             end
              if obj.RequiresFutilityThreshold
                   if ~isfield(strategyParams, 'VC_not') || ~isvector(strategyParams.VC_not) || isempty(strategyParams.VC_not)
                       error('%s: Futility enabled (RequiresFutilityThreshold=true), but strategyParams.VC_not vector is missing or empty.', obj.Name);
                   end
              end

             MM = analysisInfo.MM;
             if iscolumn(MM); MM = MM'; end % Ensure row vector
             alpha = strategyParams.alpha;
             K = numel(MM);

             % --- Calculate Detection Thresholds (VC_MSC logic) ---
             VC_detection = nan(1, K);
             for i = 1:K
                 m_val = MM(i);
                 if m_val > 1
                     VC_detection(i) = max(0, min(1, 1 - alpha^(1 / (m_val - 1))));
                 else
                     VC_detection(i) = NaN;
                     warning('%s:calculateThresholds', 'Detection threshold undefined for M=1 at stage %d. Set to NaN.', obj.Name, i);
                 end
             end
             thresholds.detection = VC_detection;
             thresholds.NDC = strategyParams.NDC;

             % --- Include Futility Thresholds (if enabled) ---
             if obj.RequiresFutilityThreshold
                 VC_not_full = strategyParams.VC_not;
                 valid_MM = MM(~isnan(MM));
                 if isempty(valid_MM)
                     warning('%s:calculateThresholds', 'Cannot determine max M for VC_not check as all MM values are NaN.', obj.Name);
                     thresholds.futility = nan(1, K);
                 else
                     max_MM = max(valid_MM);
                     if numel(VC_not_full) < max_MM
                         error('%s: Provided VC_not vector length (%d) is less than the maximum valid M required (%d).', obj.Name, numel(VC_not_full), max_MM);
                     end
                     try
                         valid_MM_indices = ~isnan(MM);
                         thresholds.futility = nan(1, K);
                         valid_MM_values = MM(valid_MM_indices);
                         if any(valid_MM_values < 1) || any(valid_MM_values > numel(VC_not_full))
                            error('%s: MM contains values outside the valid index range [1, %d] for VC_not.', obj.Name, numel(VC_not_full));
                         end
                         thresholds.futility(valid_MM_indices) = VC_not_full(valid_MM_values);
                         thresholds.futility = thresholds.futility(:)';
                     catch ME_index
                          error('%s: Error indexing VC_not with MM values. Ensure MM contains valid indices for VC_not. Error: %s', obj.Name, ME_index.message);
                     end
                 end
             else
                  thresholds.futility = []; % Indicate futility is not used
             end
        end % calculateThresholds


        function [decisionSequence, stoppingStage] = runTest(obj, statSequence, thresholds)
             % Runs the ETS test (NDC or Pnv) on a sequence of test statistics.
             % Args:
             %   statSequence: Test statistics (e.g., MSC) (bins x stages x channels). Stages correspond to MM points.
             %   thresholds: Struct from calculateThresholds (.detection, .futility, .NDC).
             % Returns:
             %   decisionSequence: Decisions (1=Detect, -1=Futility, 0=Continue) (bins x stages x channels).
             %   stoppingStage: Stage number of decision (bins x channels).

             if ~isfield(thresholds, 'detection') || ~isfield(thresholds, 'NDC')
                 error('%s: Threshold struct missing detection or NDC fields.', obj.Name);
             end
             if isempty(statSequence); decisionSequence=[]; stoppingStage=[]; return; end

             [nBins, K, nChannels] = size(statSequence);
             VC_det = thresholds.detection(:)';
             NDC = thresholds.NDC;

             if numel(VC_det) ~= K; error('%s: Detection threshold length (%d) mismatch with stages (%d).', obj.Name, numel(VC_det), K); end

             VC_fut = [];
             if obj.RequiresFutilityThreshold
                  if ~isfield(thresholds, 'futility'); error('%s: Futility enabled but futility thresholds field missing.', obj.Name); end
                  VC_fut = thresholds.futility(:)';
                   if ~isempty(VC_fut) && numel(VC_fut) ~= K; error('%s: Futility threshold length (%d) mismatch with stages (%d).', obj.Name, numel(VC_fut), K); end
             end

             decisionSequence = zeros(nBins, K, nChannels);
             stoppingStage = ones(nBins, nChannels) * K;

             for chan = 1:nChannels
                 for freqBin = 1:nBins
                     consecutiveDetections = 0;
                     alreadyStopped = false;
                     for k = 1:K % Loop through stages
                         if alreadyStopped
                              if k > 1; decisionSequence(freqBin, k, chan) = decisionSequence(freqBin, k-1, chan); end
                              continue;
                         end

                         currentStat = statSequence(freqBin, k, chan);

                         if isnan(currentStat)
                             decisionSequence(freqBin, k:K, chan) = NaN;
                             stoppingStage(freqBin, chan) = k; alreadyStopped = true; continue;
                         end

                         % Check Detection
                         isDetect = false;
                         if ~isnan(VC_det(k)) && (currentStat > VC_det(k)); isDetect = true; end

                         if isDetect; consecutiveDetections = consecutiveDetections + 1;
                         else; consecutiveDetections = 0; end

                         if consecutiveDetections >= NDC
                             decisionSequence(freqBin, k:K, chan) = 1;
                             stoppingStage(freqBin, chan) = k; alreadyStopped = true; continue;
                         end

                         % Check Futility
                         if obj.RequiresFutilityThreshold && ~isempty(VC_fut) && ~isnan(VC_fut(k))
                             isFutile = (currentStat < VC_fut(k));
                             if isFutile
                                 decisionSequence(freqBin, k:K, chan) = -1;
                                 stoppingStage(freqBin, chan) = k; alreadyStopped = true; continue;
                             end
                         end
                     end % stage loop (k)
                 end % frequency loop (freqBin)
             end % channel loop (chan)
        end % runTest method

     end % methods
end % classdef