function aggregatedResultsMap = runExperimentBatch(config, dtl, ppc, ordc, tester)
% Runs the experiment loop for given configuration and objects.
% Aggregates results per unique combination of ORD parameters across subjects/stimuli.

fprintf('[%s] Starting experiment batch processing...\n', datetime);

aggregatedResultsMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
totalIterations = numel(config.subjects) * numel(config.stimuli) * numel(config.ord_startWindows) * numel(config.ord_K_stages);
processedCount = 0;

startTime = tic; % Start timer

for iSub = 1:numel(config.subjects)
    subj = config.subjects(iSub);
    for iStim = 1:numel(config.stimuli)
        stim = config.stimuli(iStim);

        % 1. Load Data
        success = dtl.loadExperiment(subj, stim);
        if ~success, warning('Skipping Subj %d, Stim %d: Load error.', subj, stim); dtl.clearData(); continue; end
        rawSignals = dtl.getRawSignals();
        currentMaxWindow = dtl.currentDuration;
        currentFs = dtl.fs;
        currentNChannels = dtl.nChannels;

        if ppc.fs ~= currentFs, ppc = PreProcessor('fs', currentFs); end
        if ordc.fs ~= currentFs, ordc = ORDCalculator(currentFs); end

        % 2. Preprocess Data
        processedSignals = rawSignals;
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
                 processedCount = processedCount + 1; % Increment counter

                 % --- ETA Calculation and Progress ---
                 if mod(processedCount, 50) == 0 || processedCount == totalIterations
                     elapsedTime = toc(startTime);
                     progressStr = sprintf('  Processed %d / %d iterations...', processedCount, totalIterations);
                     if elapsedTime > 60 && processedCount > 0
                         timePerIter = elapsedTime / processedCount;
                         remainingIter = totalIterations - processedCount;
                         etaSeconds = remainingIter * timePerIter;
                         etaMinutes = floor(etaSeconds / 60);
                         etaSecRem = round(mod(etaSeconds, 60));
                         etaStr = sprintf(' ETA: %d min %02d sec', etaMinutes, etaSecRem);
                         progressStr = [progressStr, etaStr]; %#ok<AGROW>
                     end
                     fprintf('%s\n', progressStr);
                 end
                 % --- End ETA ---

                 if startWin >= currentMaxWindow, continue; end

                 % 4. Calculate MSC
                 [msc, epochs, M, K_actual] = ordc.calculateMSC(processedSignals, ...
                                                                 'method', config.ord_method, ...
                                                                 'startWindow', startWin, ...
                                                                 'K', K_val, ...
                                                                 'maxWindow', currentMaxWindow);
                 if isempty(msc) || K_actual < 1 || M < 2, continue; end

                 paramKey = sprintf('SW%d_K%d_Mthd%s', startWin, K_val, config.ord_method);
                 currentOrdParams = struct('startWindow', startWin, 'K', K_val, 'method', config.ord_method, 'M', M);
                 % Add subject/stimulus info to ORD params for potential use later if needed
                 currentOrdParams.subject = subj; % Store subj/stim associated with this data point before aggregation
                 currentOrdParams.stimulus = stim;

                 % 5. Run Tester
                 [decisionResults, ~, ~] = tester.runBetaCGST(msc, M, K_actual);

                 % 6. Extract counts
                 tp_matrix = decisionResults.TP(1:dtl.nSignalFrequencies, :, :);
                 fn_matrix = decisionResults.FN(1:dtl.nSignalFrequencies, :, :);
                 fp_matrix = decisionResults.FP(dtl.nSignalFrequencies+1:end, :, :);
                 tn_matrix = decisionResults.TN(dtl.nSignalFrequencies+1:end, :, :);
                 tp_this = sum(any(tp_matrix, 2), [1, 3]); fn_this = sum(any(fn_matrix, 2), [1, 3]);
                 fp_this = sum(any(fp_matrix, 2), [1, 3]); tn_this = sum(any(tn_matrix, 2), [1, 3]);
                 nSigOpps_this = dtl.nSignalFrequencies * currentNChannels;
                 nNoiseOpps_this = dtl.nNoiseFrequencies * currentNChannels;

                 % 7. Aggregate into the map (SUMMING across subjects/stimuli)
                 if aggregatedResultsMap.isKey(paramKey)
                     aggData = aggregatedResultsMap(paramKey);
                     aggData.sumTP = aggData.sumTP + tp_this; aggData.sumFN = aggData.sumFN + fn_this;
                     aggData.sumFP = aggData.sumFP + fp_this; aggData.sumTN = aggData.sumTN + tn_this;
                     aggData.totalSignalOpps = aggData.totalSignalOpps + nSigOpps_this;
                     aggData.totalNoiseOpps = aggData.totalNoiseOpps + nNoiseOpps_this;
                     % Keep track of M values encountered for this ORD setting
                     if ~isfield(aggData, 'MList'), aggData.MList = []; end
                     aggData.MList = [aggData.MList, M];
                     aggregatedResultsMap(paramKey) = aggData;
                 else
                      aggregatedResultsMap(paramKey) = struct(...
                          'ordParams', rmfield(currentOrdParams, {'subject','stimulus'}),... % Store only ORD params + M
                          'sumTP', tp_this, 'sumFN', fn_this, 'sumFP', fp_this, 'sumTN', tn_this,...
                          'totalSignalOpps', nSigOpps_this, 'totalNoiseOpps', nNoiseOpps_this, ...
                          'MList', M ... % Initialize M list
                          );
                 end
                 clear msc decisionResults epochs tp_matrix fn_matrix fp_matrix tn_matrix;
             end % K loop
        end % startWindow loop
        dtl.clearData();
        clear rawSignals processedSignals;
    end % Stimulus loop
end % Subject loop

% Final progress update
if mod(processedCount, 50) ~= 0 && processedCount == totalIterations
     fprintf('  Processed %d / %d iterations...\n', processedCount, totalIterations);
end
fprintf('[%s] Experiment batch processing finished. Total time: %.2f seconds.\n', datetime, toc(startTime));

end