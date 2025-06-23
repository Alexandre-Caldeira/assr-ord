% --- FILE START: SequentialTester.m ---
classdef SequentialTester < handle
    % Orchestrates the execution of sequential tests using a specified strategy.
    % Manages test parameters, runs tests on bulk data, and calculates performance metrics.

    properties
        Strategy % Handle to the TestStrategy object to be used. <<< TYPE HINT REMOVED
        DataManager DataManager % Handle to DataManager for metadata (Fs, frequencies, mode).
        MSCAnalyzer MSCAnalyzer % Handle to MSCAnalyzer for accessing MSC data and epoch info.

        % Results Storage
        TestDecisions % Cell array {stim/mean, subj/std, param1, param2} -> decisions(bin, stage, channel)
        StoppingStages % Cell array {stim/mean, subj/std, param1, param2} -> stoppingStage(bin, channel)
        StrategyParamsArray % Cell array storing the parameters used for each test run.
        PerformanceMetrics % Struct storing performance metrics (e.g., TPR, FPR). Calculated by calculatePerformance.
    end

    methods
        function obj = SequentialTester(strategy, dataManager, mscAnalyzer)
            % Constructor: Initializes the tester with necessary components.
            arguments
                strategy TestStrategy % <<< Keep type hint for constructor input validation
                dataManager DataManager
                mscAnalyzer MSCAnalyzer
            end
            % Assign the concrete strategy object passed in
            obj.Strategy = strategy;
            obj.DataManager = dataManager;
            obj.MSCAnalyzer = mscAnalyzer;
            fprintf('SequentialTester: Initialized with Strategy: %s.\n', obj.Strategy.Name); % Can now safely access strategy name
        end

        function obj = runBulkTest(obj, strategyParamsInput)
            % Runs the configured test strategy on all exams and parameter sets
            % defined in the MSCAnalyzer.
            % Args:
            %   strategyParamsInput: Either a scalar struct (if params are the same for all tests)
            %                        or a cell array matching the dimensions of MSCAnalyzer.mscValues
            %                        ({stim/mean, subj/std, p1, p2}), where each cell contains the
            %                        struct needed by the Strategy's calculateThresholds method
            %                        (e.g., struct('alpha', 0.05, 'NDC', 3)).

            bulkMSCData = obj.MSCAnalyzer.getBulkMSCData();
            epochDefs = obj.MSCAnalyzer.epochDefinitions; % Get corresponding epoch info
            winPerStage = obj.MSCAnalyzer.windowsPerStage; % Get M info

            if isempty(bulkMSCData)
                error('SequentialTester: MSC data is empty in MSCAnalyzer. Run MSC analysis first.');
            end

            [nDim1, nDim2, nParam1, nParam2] = size(bulkMSCData);
            dataSize = [nDim1, nDim2, nParam1, nParam2];

            % --- Validate and Store Strategy Parameters ---
            useScalarParams = false;
            if isstruct(strategyParamsInput) && isscalar(strategyParamsInput)
                obj.StrategyParamsArray = repmat({strategyParamsInput}, dataSize); % Expand scalar to cell array
                useScalarParams = true;
                fprintf('SequentialTester: Using scalar strategy parameters for all tests.\n');
            elseif iscell(strategyParamsInput) && isequal(size(strategyParamsInput), dataSize)
                obj.StrategyParamsArray = strategyParamsInput;
                fprintf('SequentialTester: Using provided cell array of strategy parameters.\n');
            else
                error('SequentialTester: strategyParamsInput must be a scalar struct or a cell array matching MSC data dimensions [%d x %d x %d x %d].', dataSize);
            end

            % --- Preallocate Results ---
            obj.TestDecisions = cell(dataSize);
            obj.StoppingStages = cell(dataSize);

            fprintf('SequentialTester: Starting bulk test execution (%s Strategy)...\n', obj.Strategy.Name);
            totalTests = numel(bulkMSCData);
            testCount = 0;

            % --- Loop Through Data and Parameter Sets ---
            for idx1 = 1:nDim1 % Stim/Mean
                for idx2 = 1:nDim2 % Subj/Std
                    for p1Idx = 1:nParam1
                        for p2Idx = 1:nParam2
                            testCount = testCount + 1;
                            fprintf('  Running Test %d/%d (Idx1:%d, Idx2:%d, P1:%d, P2:%d)...', ...
                                    testCount, totalTests, idx1, idx2, p1Idx, p2Idx);

                            mscData = bulkMSCData{idx1, idx2, p1Idx, p2Idx};
                            epochs = epochDefs{idx1, idx2, p1Idx, p2Idx};
                            M_vals = winPerStage{idx1, idx2, p1Idx, p2Idx};
                            currentStrategyParams = obj.StrategyParamsArray{idx1, idx2, p1Idx, p2Idx};

                            % --- Validate Inputs for this Test ---
                            if isempty(mscData) || isempty(epochs) || isempty(M_vals) || isempty(currentStrategyParams)
                                fprintf(' Skipping (missing data/params).\n');
                                obj.TestDecisions{idx1, idx2, p1Idx, p2Idx} = []; % Mark as skipped
                                obj.StoppingStages{idx1, idx2, p1Idx, p2Idx} = [];
                                continue;
                            end
                            if numel(epochs) < 2
                                 fprintf(' Skipping (invalid epoch definition).\n');
                                 continue;
                            end

                            % --- Prepare Analysis Info for Strategy ---
                            analysisInfo.K = size(mscData, 2); % Number of stages from MSC data
                            if analysisInfo.K == 0
                                 fprintf(' Skipping (MSC data has 0 stages).\n');
                                 continue;
                            end
                            analysisInfo.M = M_vals; % Windows per stage (scalar or vector)
                            % Add MM for ETS strategies (window count at each test point)
                             if numel(epochs) > 1
                                analysisInfo.MM = epochs(2:end) - 1; % M value at the end of stage k
                             else
                                 analysisInfo.MM = []; % Or handle error appropriately
                                  fprintf(' Skipping (cannot determine MM from epochs).\n');
                                  continue;
                             end


                             % --- Calculate Thresholds ---
                            try
                                thresholds = obj.Strategy.calculateThresholds(currentStrategyParams, analysisInfo);
                                fprintf(' Thresh OK...');
                            catch ME
                                warning('SequentialTester:runBulkTest', '\n    Error calculating thresholds for Idx1:%d, Idx2:%d, P1:%d, P2:%d: %s. Skipping test.', idx1, idx2, p1Idx, p2Idx, ME.message);
                                 obj.TestDecisions{idx1, idx2, p1Idx, p2Idx} = [];
                                 obj.StoppingStages{idx1, idx2, p1Idx, p2Idx} = [];
                                continue; % Skip to next test
                            end

                             % --- Run the Test ---
                             try
                                [decisions, stopping] = obj.Strategy.runTest(mscData, thresholds);
                                obj.TestDecisions{idx1, idx2, p1Idx, p2Idx} = decisions;
                                obj.StoppingStages{idx1, idx2, p1Idx, p2Idx} = stopping;
                                fprintf(' Test OK.\n');
                             catch ME_run
                                  warning('SequentialTester:runBulkTest','\n    Error running test for Idx1:%d, Idx2:%d, P1:%d, P2:%d: %s. Skipping test.', idx1, idx2, p1Idx, p2Idx, ME_run.message);
                                   obj.TestDecisions{idx1, idx2, p1Idx, p2Idx} = [];
                                   obj.StoppingStages{idx1, idx2, p1Idx, p2Idx} = [];
                                   continue; % Skip to next test
                             end

                        end % p2Idx
                    end % p1Idx
                end % idx2
            end % idx1
            fprintf('SequentialTester: Bulk testing complete.\n');
        end % runBulkTest

        function obj = calculatePerformance(obj)
            % Calculates performance metrics (TPR, FPR, TNR, FNR) based on stored test decisions.
            % Assumes Simulation mode where ground truth (signal/noise frequencies) is known.

            if isempty(obj.TestDecisions)
                 warning('SequentialTester: Cannot calculate performance. Run runBulkTest first.');
                 obj.PerformanceMetrics = [];
                 return;
            end
            if ~strcmp(obj.DataManager.Mode, 'sim')
                 warning('SequentialTester: Performance calculation requires simulation mode for ground truth. Results might be inaccurate.');
                 % Proceed but metrics might be meaningless without known ground truth
            end

            fprintf('SequentialTester: Calculating performance metrics...\n');
            decisionsCell = obj.TestDecisions;
            stoppingCell = obj.StoppingStages;
            if isempty(decisionsCell) % Check if TestDecisions is populated
                warning('SequentialTester: TestDecisions cell array is empty. No performance calculated.');
                obj.PerformanceMetrics = [];
                return;
            end
            [nDim1, nDim2, nParam1, nParam2] = size(decisionsCell); % Dims match noise means/stds/param combos

            % --- Get Ground Truth Frequencies ---
            fs = obj.DataManager.Fs;
            nfft = obj.DataManager.getNFFT();
            if nfft == 0; error('SequentialTester: NFFT is zero, cannot map frequencies to bins.'); end
            freqResolution = fs / nfft;

            sigFreqs_Hz = obj.DataManager.SimSignalFrequencies;
            noiseFreqs_Hz = obj.DataManager.SimNoiseFrequencies;

            % Convert Hz to approximate bin indices (1-based)
            sigFreqBins = round(sigFreqs_Hz / freqResolution) + 1;
            noiseFreqBins = round(noiseFreqs_Hz / freqResolution) + 1;

            % Ensure bins are within valid range
            maxBin = obj.DataManager.getNumBins();
            sigFreqBins = sigFreqBins(sigFreqBins >= 1 & sigFreqBins <= maxBin);
            noiseFreqBins = noiseFreqBins(noiseFreqBins >= 1 & noiseFreqBins <= maxBin);
            sigFreqBins = unique(sigFreqBins);
            noiseFreqBins = unique(noiseFreqBins);
            noiseFreqBins = setdiff(noiseFreqBins, sigFreqBins); % Remove overlap

            nSignalFreqs = numel(sigFreqBins);
            nNoiseFreqs = numel(noiseFreqBins);

            if nSignalFreqs == 0; warning('SequentialTester: No valid signal frequency bins found after mapping/validation.'); end
            if nNoiseFreqs == 0; warning('SequentialTester: No valid noise frequency bins found after mapping/validation.'); end
            fprintf('  Using %d signal bins, %d noise bins.\n', nSignalFreqs, nNoiseFreqs);

            % --- Store results per parameter set ---
            paramResults = repmat(struct('TP', 0, 'FP', 0, 'TN', 0, 'FN', 0, ...
                                         'TotalSignalOpportunities', 0, 'TotalNoiseOpportunities', 0), ...
                                  [nParam1, nParam2]);

            % --- Aggregate results ---
            for p1Idx = 1:nParam1
                for p2Idx = 1:nParam2
                    fprintf('  Aggregating performance for ParamSet (P1:%d, P2:%d)...\n', p1Idx, p2Idx);
                    tp_count = 0; fp_count = 0; tn_count = 0; fn_count = 0;
                    total_sig_opp = 0; total_noise_opp = 0;

                    for idx1 = 1:nDim1 % Noise Mean index
                        for idx2 = 1:nDim2 % Noise Std index
                            decisions = decisionsCell{idx1, idx2, p1Idx, p2Idx};
                            stopping = stoppingCell{idx1, idx2, p1Idx, p2Idx};

                            if isempty(decisions) || isempty(stopping); continue; end

                            [nBins, K, nChannels] = size(decisions);
                            if isempty(nChannels) || nChannels == 0 || K == 0; continue; end

                            % Get final decision state
                            final_decision = zeros(nBins, nChannels);
                            for chan = 1:nChannels
                                for bin = 1:nBins
                                     % Check if stopping indices are valid for this bin/chan
                                     if bin > size(stopping, 1) || chan > size(stopping, 2)
                                         % warning('SequentialTester: Stopping stage index (%d, %d) out of bounds for size [%d x %d]. Using K as stop stage.', bin, chan, size(stopping, 1), size(stopping, 2));
                                         stop_k = K; % Default to last stage if index invalid
                                     else
                                         stop_k = stopping(bin, chan);
                                     end

                                    % Handle cases where stopping might be NaN or invalid
                                    if isnan(stop_k) || stop_k < 1 || stop_k > K
                                        stop_k = K; % Default to checking last stage if stop info invalid
                                    else
                                        stop_k = round(min(stop_k, K)); % Ensure integer within bounds [1, K]
                                    end

                                    % Check if decision indices are valid
                                     if bin > size(decisions, 1) || stop_k > size(decisions, 2) || chan > size(decisions, 3)
                                          warning('SequentialTester: Decision index (%d, %d, %d) out of bounds for size [%d x %d x %d]. Skipping bin/channel.', bin, stop_k, chan, size(decisions, 1), size(decisions, 2), size(decisions, 3));
                                         final_decision(bin, chan) = NaN; % Mark as invalid
                                         continue;
                                     end

                                    final_decision_val = decisions(bin, stop_k, chan);
                                    % Treat NaN decision as "no detection" (FN for signal, TN for noise)
                                    if isnan(final_decision_val)
                                        final_decision(bin, chan) = 0; % Or -1 if preferred interpretation
                                    else
                                        final_decision(bin, chan) = final_decision_val;
                                    end
                                end % bin loop
                            end % chan loop (for final decision)


                            % Count outcomes vs ground truth
                            for chan = 1:nChannels
                                % Signal Frequencies
                                if nSignalFreqs > 0 && ~isempty(sigFreqBins)
                                    valid_sig_bins = sigFreqBins(sigFreqBins <= nBins); % Ensure bins exist in decisions
                                    if ~isempty(valid_sig_bins)
                                         chan_decisions_sig = final_decision(valid_sig_bins, chan);
                                         tp_count = tp_count + sum(chan_decisions_sig == 1);
                                         fn_count = fn_count + sum(chan_decisions_sig ~= 1); % Includes 0, -1, NaN
                                         total_sig_opp = total_sig_opp + numel(valid_sig_bins);
                                    end
                                end
                                % Noise Frequencies
                                  if nNoiseFreqs > 0 && ~isempty(noiseFreqBins)
                                     valid_noise_bins = noiseFreqBins(noiseFreqBins <= nBins); % Ensure bins exist
                                     if ~isempty(valid_noise_bins)
                                         chan_decisions_noise = final_decision(valid_noise_bins, chan);
                                         fp_count = fp_count + sum(chan_decisions_noise == 1);
                                         tn_count = tn_count + sum(chan_decisions_noise ~= 1); % Includes 0, -1, NaN
                                         total_noise_opp = total_noise_opp + numel(valid_noise_bins);
                                     end
                                end
                            end % channel loop (for counting)
                        end % idx2 loop (Std)
                    end % idx1 loop (Mean)

                    % Store aggregated counts
                    paramResults(p1Idx, p2Idx).TP = tp_count;
                    paramResults(p1Idx, p2Idx).FP = fp_count;
                    paramResults(p1Idx, p2Idx).TN = tn_count;
                    paramResults(p1Idx, p2Idx).FN = fn_count;
                    paramResults(p1Idx, p2Idx).TotalSignalOpportunities = total_sig_opp;
                    paramResults(p1Idx, p2Idx).TotalNoiseOpportunities = total_noise_opp;
                end % p2Idx loop
            end % p1Idx loop

            obj.PerformanceMetrics.PerParameter = paramResults;

            % --- Calculate Overall Rates ---
             try % Wrap in try-catch for safety
                 all_TP = sum([paramResults.TP]); all_FP = sum([paramResults.FP]);
                 all_TN = sum([paramResults.TN]); all_FN = sum([paramResults.FN]);
                 all_TotalSig = sum([paramResults.TotalSignalOpportunities]);
                 all_TotalNoise = sum([paramResults.TotalNoiseOpportunities]);

                 obj.PerformanceMetrics.Overall = struct(); % Initialize struct

                 if all_TotalSig > 0
                     obj.PerformanceMetrics.Overall.TPR = 100 * all_TP / all_TotalSig;
                     obj.PerformanceMetrics.Overall.FNR = 100 * all_FN / all_TotalSig;
                 else
                     obj.PerformanceMetrics.Overall.TPR = NaN; obj.PerformanceMetrics.Overall.FNR = NaN;
                 end
                 if all_TotalNoise > 0
                     obj.PerformanceMetrics.Overall.FPR = 100 * all_FP / all_TotalNoise;
                     obj.PerformanceMetrics.Overall.TNR = 100 * all_TN / all_TotalNoise;
                 else
                      obj.PerformanceMetrics.Overall.FPR = NaN; obj.PerformanceMetrics.Overall.TNR = NaN;
                 end
                 if (all_TotalSig + all_TotalNoise) > 0
                    obj.PerformanceMetrics.Overall.Accuracy = 100 * (all_TP + all_TN) / (all_TotalSig + all_TotalNoise);
                 else
                     obj.PerformanceMetrics.Overall.Accuracy = NaN;
                 end

                 fprintf('SequentialTester: Performance calculation complete.\n');
                 disp('Overall Performance Metrics:');
                 disp(obj.PerformanceMetrics.Overall);

            catch ME_perf
                warning('SequentialTester:calculatePerformance', 'Error calculating overall performance metrics: %s', ME_perf.message);
                obj.PerformanceMetrics.Overall = struct('Error', ME_perf.message);
            end

         end % calculatePerformance method

         function obj = calculateAverageStoppingTime(obj, p)
             % Calculates average stopping time (in stages or seconds) per parameter set.
             arguments
                 obj
                 p.Units {mustBeMember(p.Units, {'stages', 'seconds'})} = 'stages'
             end

             avgStopTime = []; % Initialize output
             if isempty(obj.StoppingStages)
                 warning('SequentialTester: Cannot calculate stopping time. Run runBulkTest first.');
                 return;
             end

             stoppingCell = obj.StoppingStages;
             [nDim1, nDim2, nParam1, nParam2] = size(stoppingCell);
             avgStopTime = nan(nParam1, nParam2); % Store avg stop time per param set

             fprintf('SequentialTester: Calculating average stopping time (Units: %s)...\n', p.Units);

             for p1Idx = 1:nParam1
                 for p2Idx = 1:nParam2
                     allStopsForParamSet = []; % Collect stopping stages for this param set

                     for idx1 = 1:nDim1
                         for idx2 = 1:nDim2
                             stopping = stoppingCell{idx1, idx2, p1Idx, p2Idx};
                             if ~isempty(stopping)
                                 allStopsForParamSet = [allStopsForParamSet; stopping(:)]; % Flatten and append
                             end
                         end
                     end

                     if ~isempty(allStopsForParamSet)
                         validStops = allStopsForParamSet(~isnan(allStopsForParamSet)); % Exclude NaNs
                         if ~isempty(validStops)
                             if strcmp(p.Units, 'stages')
                                 avgStopTime(p1Idx, p2Idx) = mean(validStops);
                             else % seconds
                                 warning('SequentialTester:calculateAverageStoppingTime', 'Calculating stopping time in seconds requires robust mapping of stages to time, which is complex in bulk mode. Returning average stage number instead.');
                                 avgStopTime(p1Idx, p2Idx) = mean(validStops);
                                 % TODO: Implement seconds calculation if epoch durations can be reliably retrieved for *each* test instance.
                             end
                         end % ~isempty(validStops)
                     end % ~isempty(allStopsForParamSet)
                 end % p2Idx
             end % p1Idx
              fprintf('SequentialTester: Average stopping time calculation finished.\n');
              % Store result within PerformanceMetrics if it exists
              if isfield(obj.PerformanceMetrics,'PerParameter') || isfield(obj.PerformanceMetrics,'Overall')
                    obj.PerformanceMetrics.AverageStoppingTime = avgStopTime;
              else
                   obj.PerformanceMetrics = struct('AverageStoppingTime',avgStopTime);
              end
         end % calculateAverageStoppingTime method

    end % methods
end % classdef