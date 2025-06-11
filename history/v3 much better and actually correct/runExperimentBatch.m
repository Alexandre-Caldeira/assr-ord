function aggregatedResultsMap = runExperimentBatch(config, dtl, ppc, ordc, tester)
% Runs the experiment loop for given configuration and objects.
% Returns a map where keys are ORD param strings and values are structs
% containing SUMMED counts/times/opportunities across subjects/stimuli.

fprintf('[%s] Starting experiment batch processing...\n', datetime);

aggregatedResultsMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
paramIdx = 0;
totalIterations = numel(config.subjects) * numel(config.stimuli) * numel(config.ord_startWindows) * numel(config.ord_K_stages);
processedCount = 0;

for iSub = 1:numel(config.subjects)
    subj = config.subjects(iSub);
    for iStim = 1:numel(config.stimuli)
        stim = config.stimuli(iStim);

        % fprintf('\n[%s] Subject %d, Stimulus %d...\n', datetime, subj, stim); % Verbose

        % 1. Load Data
        success = dtl.loadExperiment(subj, stim);
        if ~success, warning('Skipping Subj %d, Stim %d: Load error.', subj, stim); dtl.clearData(); continue; end
        rawSignals = dtl.getRawSignals();
        currentMaxWindow = dtl.currentDuration;
        currentFs = dtl.fs;
        currentNChannels = dtl.nChannels; % N channels for this specific file

        % Update dependent objects if Fs changed from previous file
        if ppc.fs ~= currentFs, ppc = PreProcessor('fs', currentFs); end % Re-init with file's Fs
        if ordc.fs ~= currentFs, ordc = ORDCalculator(currentFs); end

        % 2. Preprocess Data
        processedSignals = rawSignals; % Start with raw
        if config.ppc_useZanoteli
            maxWinPPC = min(dtl.zanoteliSuggestedMMax(stim), size(processedSignals, 2));
            processedSignals = ppc.applyZanoteliPreprocessing(processedSignals, maxWinPPC);
        end
        if config.ppc_useFilter && ~isempty(processedSignals)
            processedSignals = ppc.applyAntunesFiltering(processedSignals);
        end

        if isempty(processedSignals) || size(processedSignals, 2) < 2
            warning('Skipping Subj %d, Stim %d: No data after PPC.', subj, stim);
            dtl.clearData(); clear rawSignals processedSignals; continue;
        end
        currentMaxWindow = size(processedSignals, 2);

        % 3. Inner loop for ORD / Tester parameters
        for iStart = 1:numel(config.ord_startWindows)
             startWin = config.ord_startWindows(iStart);
             for iK = 1:numel(config.ord_K_stages)
                 K_val = config.ord_K_stages(iK);
                 paramIdx = paramIdx + 1;
                 processedCount = processedCount + 1;
                 if mod(processedCount, 50) == 0 % Print progress marker every 50 iterations
                     fprintf('  Processed %d / %d iterations...\n', processedCount, totalIterations);
                 end

                 paramKey = sprintf('SW%d_K%d_Mthd%s', startWin, K_val, config.ord_method);
                 currentOrdParams = struct('startWindow', startWin, 'K', K_val, 'method', config.ord_method); % For results table later

                 if startWin >= currentMaxWindow, continue; end % Skip if start window invalid

                 % 4. Calculate MSC
                 [msc, epochs, M, K_actual] = ordc.calculateMSC(processedSignals, ...
                                                                 'method', config.ord_method, ...
                                                                 'startWindow', startWin, ...
                                                                 'K', K_val, ...
                                                                 'maxWindow', currentMaxWindow);

                 if isempty(msc) || K_actual < 1 || M < 2, continue; end % Skip if MSC failed

                 % 5. Run Tester
                 [decisionResults, ~, ~] = tester.runBetaCGST(msc, M, K_actual);

                 % 6. Extract counts for THIS subject/stim/ORD_param config
                 tp_matrix = decisionResults.TP(1:dtl.nSignalFrequencies, :, :);
                 fn_matrix = decisionResults.FN(1:dtl.nSignalFrequencies, :, :);
                 fp_matrix = decisionResults.FP(dtl.nSignalFrequencies+1:end, :, :);
                 tn_matrix = decisionResults.TN(dtl.nSignalFrequencies+1:end, :, :);

                 tp_this = sum(any(tp_matrix, 2), [1, 3]);
                 fn_this = sum(any(fn_matrix, 2), [1, 3]);
                 fp_this = sum(any(fp_matrix, 2), [1, 3]);
                 tn_this = sum(any(tn_matrix, 2), [1, 3]);

                 nSigOpps_this = dtl.nSignalFrequencies * currentNChannels;
                 nNoiseOpps_this = dtl.nNoiseFrequencies * currentNChannels;

                 % Time aggregation (can be added similarly if needed)
                 % [min_t, avg_t, n_t] = calculateDetectionTimes(tp_matrix, epochs);


                 % 7. Aggregate into the map
                 if aggregatedResultsMap.isKey(paramKey)
                     % Update existing entry
                     aggData = aggregatedResultsMap(paramKey);
                     aggData.sumTP = aggData.sumTP + tp_this;
                     aggData.sumFN = aggData.sumFN + fn_this;
                     aggData.sumFP = aggData.sumFP + fp_this;
                     aggData.sumTN = aggData.sumTN + tn_this;
                     aggData.totalSignalOpps = aggData.totalSignalOpps + nSigOpps_this;
                     aggData.totalNoiseOpps = aggData.totalNoiseOpps + nNoiseOpps_this;
                     % Update time sums/counts here
                     % aggData.sumMinDetTime = aggData.sumMinDetTime + min_t; % Need careful handling of NaN/min
                     % aggData.nMinDetTime = aggData.nMinDetTime + (1 if ~isnan(min_t));
                     % aggData.sumAvgDetTime = aggData.sumAvgDetTime + avg_t * n_t; % Sum of times
                     % aggData.nAvgDetTime = aggData.nAvgDetTime + n_t; % Total detections contributing to avg
                     aggregatedResultsMap(paramKey) = aggData;
                 else
                     % Create new entry
                      aggregatedResultsMap(paramKey) = struct(...
                          'ordParams', currentOrdParams,... % Store ORD params once
                          'sumTP', tp_this, 'sumFN', fn_this, 'sumFP', fp_this, 'sumTN', tn_this,...
                          'totalSignalOpps', nSigOpps_this, 'totalNoiseOpps', nNoiseOpps_this ...
                          % Initialize time fields here if aggregating time
                          );
                 end

                 clear msc decisionResults epochs tp_matrix fn_matrix fp_matrix tn_matrix;
             end % K loop
        end % startWindow loop

        dtl.clearData();
        clear rawSignals processedSignals;
        % fprintf('[%s] Finished Subject %d, Stimulus %d.\n', datetime, subj, stim); % Verbose
    end % Stimulus loop
end % Subject loop

fprintf('[%s] Experiment batch processing finished. Processed %d iterations.\n', datetime, processedCount);

end % End function runExperimentBatch