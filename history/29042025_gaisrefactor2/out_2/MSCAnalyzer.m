% --- FILE START: MSCAnalyzer.m ---
classdef MSCAnalyzer
    % Calculates Magnitude-Squared Coherence (MSC) based on different
    % epoching strategies using frequency-domain data.

    properties
        EpochMethod {mustBeMember(EpochMethod, {'fixedSize', 'fixedK'})} = 'fixedSize' % Method to define stages: 'fixedSize' (constant step) or 'fixedK' (fixed number of stages).
        % --- Parameters for 'fixedSize' method ---
        StartWindow (1,1) double {mustBeInteger, mustBePositive} = 1 % Starting window index for the first stage.
        WindowStepSize (1,1) double {mustBeInteger, mustBePositive} = 24 % Number of windows per stage (step size).
        % --- Parameters for 'fixedK' method ---
        NumStages (1,1) double {mustBeInteger, mustBePositive} = 5 % Desired number of stages (K).

        % --- Output Storage ---
        mscValues % Cell array {stim/mean, subj/std, param1_idx, param2_idx} -> msc(bin, stage, channel)
        epochDefinitions % Cell array {stim/mean, subj/std, param1_idx, param2_idx} -> [epoch_boundary_indices]
        windowsPerStage % Cell array {stim/mean, subj/std, param1_idx, param2_idx} -> M (scalar or vector)

        % --- Internal Parameter Storage ---
        analyzedParam1Values = []
        analyzedParam2Values = []
        analyzedParam1Name = ''
        analyzedParam2Name = ''
    end

    methods
        function obj = MSCAnalyzer(epochMethod, options)
            % Constructor: Sets epoching method and associated parameters using Name-Value pairs.
            % Args:
            %   epochMethod: 'fixedSize' or 'fixedK'.
            %   options: Name-value pairs for ('StartWindow', 'WindowStepSize', 'NumStages').
            arguments % <<< Use arguments block for validation and Name-Value pairs
                epochMethod {mustBeMember(epochMethod, {'fixedSize', 'fixedK'})}
                % Define Name-Value arguments within the options struct
                options.StartWindow (1,1) double {mustBeInteger, mustBePositive} = 1 % Default value if not provided
                options.WindowStepSize (1,1) double {mustBeInteger, mustBePositive} = 24 % Default value
                options.NumStages (1,1) double {mustBeInteger, mustBePositive} = 5 % Default value
            end % <<< End arguments block

            obj.EpochMethod = epochMethod;

            % Assign values from the 'options' struct created by the arguments block
            obj.StartWindow = options.StartWindow;
            obj.WindowStepSize = options.WindowStepSize;
            obj.NumStages = options.NumStages;

            fprintf('MSCAnalyzer: Initialized with EpochMethod=%s.\n', obj.EpochMethod);
             if strcmp(obj.EpochMethod, 'fixedSize')
                 fprintf('  Default/Specified fixedSize params: StartWindow=%d, WindowStepSize=%d\n', obj.StartWindow, obj.WindowStepSize);
             else % fixedK
                 fprintf('  Default/Specified fixedK params: StartWindow=%d, NumStages=%d\n', obj.StartWindow, obj.NumStages);
             end
        end % Constructor


        function [mscResult, epochs, M_winsPerStage] = calculateMSC(obj, freqDomainData, currentParams)
            % Calculates MSC for a single exam's FFT data using specified parameters.
            % Args:
            %   freqDomainData: FFT data (bins x windows x channels).
            %   currentParams: Struct containing parameters for this specific calculation
            %                  (e.g., .StartWindow, .WindowStepSize or .NumStages).
            % Returns:
            %   mscResult: MSC values (bins x stages x channels), NaN if error.
            %   epochs: Vector of window indices defining stage boundaries.
            %   M_winsPerStage: Number of windows in each calculated stage (scalar or vector).

             mscResult = []; epochs = []; M_winsPerStage = []; % Initialize outputs

            if isempty(freqDomainData) || all(isnan(freqDomainData(:)), 'all')
                % warning('calculateMSC: Input freqDomainData is empty or all NaN. Cannot calculate MSC.'); % Less verbose
                return;
            end

            [nBins, nWindowsTotal, nChannels] = size(freqDomainData);
             if nWindowsTotal < 2
                 warning('calculateMSC: Cannot define stages with fewer than 2 windows (found %d).', nWindowsTotal);
                 return;
             end

            % --- Determine Epoch Boundaries ---
            startWin = currentParams.StartWindow; % Use params passed for this specific calculation

            switch obj.EpochMethod % Use obj.EpochMethod to decide logic path
                case 'fixedSize'
                    stepSize = currentParams.WindowStepSize;
                    if startWin > nWindowsTotal
                         warning('calculateMSC (fixedSize): StartWindow (%d) exceeds total windows (%d).', startWin, nWindowsTotal); return;
                    end
                    epoch_ends = startWin:stepSize:nWindowsTotal;
                    epochs = unique([startWin, epoch_ends + 1]);
                    if isempty(epochs) || epochs(end) <= nWindowsTotal
                         epochs = [epochs, nWindowsTotal + 1];
                    elseif epochs(end) > nWindowsTotal + 1
                         epochs(end) = nWindowsTotal + 1;
                    end
                    epochs = epochs(epochs <= nWindowsTotal + 1);
                    epochs = unique(epochs);

                    if numel(epochs) < 2
                         warning('calculateMSC (fixedSize): Epoch definition resulted in < 2 boundaries (Start=%d, Step=%d, TotalWin=%d).', startWin, stepSize, nWindowsTotal); return;
                    end
                     numStagesActual = numel(epochs) - 1;
                     M_winsPerStage = repmat(stepSize, 1, numStagesActual);
                     lastStageEnd = epochs(end)-1;
                     lastStageStart = epochs(end-1);
                     M_last = lastStageEnd - lastStageStart + 1;
                      if M_last >= 0
                         M_winsPerStage(end) = M_last;
                     else
                         warning('calculateMSC (fixedSize): Invalid M calculated for last stage (%d). Check epoch boundaries.',M_last);
                         epochs = epochs(1:end-1);
                          if numel(epochs)<2; return; end
                          numStagesActual = numel(epochs) - 1;
                          M_winsPerStage = M_winsPerStage(1:numStagesActual);
                     end

                case 'fixedK'
                    K = currentParams.NumStages;
                     if startWin > nWindowsTotal
                         warning('calculateMSC (fixedK): StartWindow (%d) exceeds total windows (%d).', startWin, nWindowsTotal); return;
                    end
                    if K < 1
                         warning('calculateMSC (fixedK): NumStages (K=%d) must be at least 1.', K); return;
                    end
                    epochs = round(linspace(startWin, nWindowsTotal + 1, K + 1));
                    epochs = unique(epochs);
                     numStagesActual = numel(epochs) - 1;
                     if numStagesActual < K
                         warning('calculateMSC (fixedK): Linspace resulted in %d stages, less than requested K=%d (Start=%d, TotalWin=%d).', numStagesActual, K, startWin, nWindowsTotal);
                          if numStagesActual < 1; return; end
                     end
                      if epochs(end) <= nWindowsTotal
                          epochs(end) = nWindowsTotal + 1;
                           epochs = unique(epochs);
                           numStagesActual = numel(epochs) - 1;
                           if numStagesActual < 1; return; end
                      end
                     M_winsPerStage = diff(epochs);

                otherwise
                    error('Unknown EpochMethod.');
            end % End epoch determination

            % --- Validate Stage Definitions ---
            if any(M_winsPerStage <= 0)
                warning('calculateMSC: Epoch definition resulted in zero or negative windows per stage. M = [%s].', num2str(M_winsPerStage));
                 mscResult = NaN(nBins, numStagesActual, nChannels); return;
             end

            % --- Calculate MSC per Stage ---
            % fprintf('  Calculating MSC for %d stages...\n', numStagesActual); % Less verbose
            mscResult = nan(nBins, numStagesActual, nChannels);

            for k = 1:numStagesActual
                idxStart = epochs(k);
                idxEnd = epochs(k+1) - 1;

                if idxStart < 1 || idxStart > nWindowsTotal || idxEnd < idxStart || idxEnd > nWindowsTotal
                     warning('  MSC Stage %d: Invalid indices [%d, %d] for total windows %d. Skipping stage.', k, idxStart, idxEnd, nWindowsTotal);
                     continue;
                end

                windowsInThisStage = idxEnd - idxStart + 1;
                 if windowsInThisStage ~= M_winsPerStage(k)
                      warning('  MSC Stage %d: M mismatch. Calculated=%d, Expected=%d. Using Calculated=%d', k, windowsInThisStage, M_winsPerStage(k), windowsInThisStage);
                      M_winsPerStage(k) = windowsInThisStage; % Correct M
                 end
                  if windowsInThisStage <= 0
                      warning('  MSC Stage %d: Calculated M=%d windows. Skipping stage.', k, windowsInThisStage);
                      continue;
                  end

                currentEpochFFT = freqDomainData(:, idxStart:idxEnd, :);

                 if any(isnan(currentEpochFFT), 'all')
                     % warning('  MSC Stage %d: NaNs detected in input FFT data. Result for this stage will be NaN.', k); % Less verbose
                     continue;
                 end

                try
                    sum_Y = sum(currentEpochFFT, 2);
                    sum_absY_sq = sum(abs(currentEpochFFT).^2, 2);
                    denominator = windowsInThisStage * sum_absY_sq;
                    zeroThreshold = eps(class(denominator));
                    valid_den_mask = denominator > zeroThreshold;
                    msc_stage = zeros(nBins, 1, nChannels) * NaN;
                    msc_stage(valid_den_mask) = abs(sum_Y(valid_den_mask)).^2 ./ denominator(valid_den_mask);
                    msc_stage = real(msc_stage);
                    msc_stage = max(0, min(1, msc_stage));
                    mscResult(:, k, :) = msc_stage;
                catch ME_msc
                     warning('  MSC Stage %d: Error during MSC calculation: %s. Result for this stage set to NaN.', k, ME_msc.message);
                     mscResult(:, k, :) = NaN;
                 end
            end % End stage loop (k)
            % fprintf('  MSC calculation finished.\n'); % Less verbose
        end % calculateMSC


        function obj = analyzeBulkMSC(obj, signalProcessor, p)
            % Calculates MSC for bulk data using ranges of parameters.
            arguments
                obj MSCAnalyzer
                signalProcessor SignalProcessor
                p.StartWindows double = obj.StartWindow
                p.WindowStepSizes double = obj.WindowStepSize
                p.NumStages double = obj.NumStages
            end

            bulkFFTData = signalProcessor.getBulkFrequencyDomainData();
            if isempty(bulkFFTData)
                error('MSCAnalyzer: Bulk FFT data is empty in SignalProcessor.');
            end
            [nDim1, nDim2] = size(bulkFFTData);

            % --- Determine Parameter Combinations ---
            if strcmp(obj.EpochMethod, 'fixedSize')
                params1_vals = p.StartWindows(:)';
                params2_vals = p.WindowStepSizes(:)';
                obj.analyzedParam1Name = 'StartWindow';
                obj.analyzedParam2Name = 'WindowStepSize';
            elseif strcmp(obj.EpochMethod, 'fixedK')
                params1_vals = p.StartWindows(:)';
                params2_vals = p.NumStages(:)';
                 obj.analyzedParam1Name = 'StartWindow';
                 obj.analyzedParam2Name = 'NumStages';
            else
                error('Unhandled EpochMethod in bulk analysis.');
            end
            nParam1 = numel(params1_vals);
            nParam2 = numel(params2_vals);
            obj.analyzedParam1Values = params1_vals;
            obj.analyzedParam2Values = params2_vals;

            % --- Preallocate Result Cell Arrays ---
            obj.mscValues = cell(nDim1, nDim2, nParam1, nParam2);
            obj.epochDefinitions = cell(nDim1, nDim2, nParam1, nParam2);
            obj.windowsPerStage = cell(nDim1, nDim2, nParam1, nParam2);

            fprintf('MSCAnalyzer: Starting bulk MSC calculation across parameters...\n');
            fprintf('  Epoch Method: %s\n', obj.EpochMethod);
            fprintf('  Parameter 1 (%s): %s\n', obj.analyzedParam1Name, mat2str(params1_vals));
            fprintf('  Parameter 2 (%s): %s\n', obj.analyzedParam2Name, mat2str(params2_vals));
            totalCalculations = nDim1 * nDim2 * nParam1 * nParam2;
            calcCount = 0;

            % --- Loop through Data and Parameters ---
            for idx1 = 1:nDim1 % Stim/Mean
                for idx2 = 1:nDim2 % Subj/Std
                    fftData = []; % Ensure defined
                    if idx1 <= size(bulkFFTData,1) && idx2 <= size(bulkFFTData,2)
                        fftData = bulkFFTData{idx1, idx2};
                    end

                    if isempty(fftData)
                         calcCount = calcCount + (nParam1 * nParam2); % Increment count even for skipped data
                         % fprintf('  Skipping Index1=%d, Index2=%d (empty FFT data). Progress: %d/%d\n', idx1, idx2, calcCount, totalCalculations); % Less verbose
                         continue;
                     end

                    for p1Idx = 1:nParam1
                        for p2Idx = 1:nParam2
                            calcCount = calcCount + 1;
                            % fprintf('  Calculating MSC %d/%d (Idx1:%d, Idx2:%d, P1:%d, P2:%d)...\n', ...
                            %         calcCount, totalCalculations, idx1, idx2, p1Idx, p2Idx); % Very verbose

                            % --- Set Parameters for this specific calculation ---
                            currentParams = struct();
                            currentParams.StartWindow = params1_vals(p1Idx);
                            if strcmp(obj.EpochMethod, 'fixedSize')
                                currentParams.WindowStepSize = params2_vals(p2Idx);
                            else % fixedK
                                currentParams.NumStages = params2_vals(p2Idx);
                            end

                            % --- Calculate MSC ---
                            [msc, epochs, M_wins] = obj.calculateMSC(fftData, currentParams);

                            % --- Store results ---
                             obj.mscValues{idx1, idx2, p1Idx, p2Idx} = msc;
                             obj.epochDefinitions{idx1, idx2, p1Idx, p2Idx} = epochs;
                             obj.windowsPerStage{idx1, idx2, p1Idx, p2Idx} = M_wins;

                        end % p2Idx
                    end % p1Idx
                     if mod(calcCount, round(totalCalculations/10))==0 || calcCount==totalCalculations % Progress update
                          fprintf('  MSC Calculation Progress: %d / %d\n', calcCount, totalCalculations);
                     end
                end % idx2
            end % idx1
             fprintf('MSCAnalyzer: Bulk MSC calculation complete.\n');
        end % analyzeBulkMSC


        % --- Access MSC Data ---
         function data = getMSCData(obj, idx1, idx2, p1Idx, p2Idx)
             % Access MSC for a specific exam and parameter set index.
             s = size(obj.mscValues);
             if isempty(s) || numel(s) ~= 4
                 warning('MSCAnalyzer: Bulk MSC values not calculated yet or have unexpected dimensions.'); data = []; return;
             end
             if idx1 >= 1 && idx1 <= s(1) && idx2 >= 1 && idx2 <= s(2) && p1Idx >= 1 && p1Idx <= s(3) && p2Idx >= 1 && p2Idx <= s(4)
                 data = obj.mscValues{idx1, idx2, p1Idx, p2Idx};
             else
                 warning('MSCAnalyzer: Indices (%d, %d, %d, %d) out of bounds for MSC data [%d x %d x %d x %d].', ...
                         idx1, idx2, p1Idx, p2Idx, s(1), s(2), s(3), s(4));
                 data = [];
             end
         end
          function bulkData = getBulkMSCData(obj)
             % Returns the entire cell array of MSC values.
             bulkData = obj.mscValues;
         end
         function epochs = getEpochs(obj, idx1, idx2, p1Idx, p2Idx)
            % Access epoch definitions for a specific exam and parameter set index.
              s = size(obj.epochDefinitions);
              if isempty(s) || numel(s) ~= 4
                 warning('MSCAnalyzer: Epoch definitions not calculated yet or have unexpected dimensions.'); epochs = []; return;
              end
             if idx1 >= 1 && idx1 <= s(1) && idx2 >= 1 && idx2 <= s(2) && p1Idx >= 1 && p1Idx <= s(3) && p2Idx >= 1 && p2Idx <= s(4)
                 epochs = obj.epochDefinitions{idx1, idx2, p1Idx, p2Idx};
             else
                 warning('MSCAnalyzer: Indices (%d, %d, %d, %d) out of bounds for epoch definitions [%d x %d x %d x %d].', ...
                          idx1, idx2, p1Idx, p2Idx, s(1), s(2), s(3), s(4));
                 epochs = [];
             end
         end
         function M = getWindowsPerStage(obj, idx1, idx2, p1Idx, p2Idx)
            % Access M (windows per stage) for a specific exam and parameter set index.
              s = size(obj.windowsPerStage);
              if isempty(s) || numel(s) ~= 4
                 warning('MSCAnalyzer: Windows per stage not calculated yet or have unexpected dimensions.'); M = []; return;
              end
             if idx1 >= 1 && idx1 <= s(1) && idx2 >= 1 && idx2 <= s(2) && p1Idx >= 1 && p1Idx <= s(3) && p2Idx >= 1 && p2Idx <= s(4)
                 M = obj.windowsPerStage{idx1, idx2, p1Idx, p2Idx};
             else
                  warning('MSCAnalyzer: Indices (%d, %d, %d, %d) out of bounds for windows per stage [%d x %d x %d x %d].', ...
                           idx1, idx2, p1Idx, p2Idx, s(1), s(2), s(3), s(4));
                 M = [];
             end
         end

    end % methods
end % classdef