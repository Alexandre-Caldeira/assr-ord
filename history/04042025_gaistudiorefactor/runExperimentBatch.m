function aggregatedResultsMap = runExperimentBatch(config, dtl, ppc, ordc, tester)
% Runs the experiment loop for given configuration and objects.
% Aggregates results per unique combination of ORD parameters across subjects/stimuli.

fprintf('[%s] Starting experiment batch processing...\n', datetime);

aggregatedResultsMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
processedCount = 0;

% --- Conditional Calculation of Total Iterations ---
if strcmp(config.ord_method, 'fixedK')
    totalIterations = numel(config.subjects) * numel(config.stimuli) * numel(config.ord_startWindows) * numel(config.ord_K_stages);
    iterParam = config.ord_K_stages;
    paramName = 'K';
elseif strcmp(config.ord_method, 'fixedStep')
    totalIterations = numel(config.subjects) * numel(config.stimuli) * numel(config.ord_startWindows) * numel(config.ord_stepSize);
    iterParam = config.ord_stepSize;
    paramName = 'stepSize';
else
    error('Invalid config.ord_method: %s', config.ord_method);
end
% ---

startTime = tic; % Start timer

for iSub = 1:numel(config.subjects)
    subj = config.subjects(iSub);
    for iStim = 1:numel(config.stimuli)
        stim = config.stimuli(iStim);

        % 1. Load Data
        success = dtl.loadExperiment(subj, stim);
        if ~success, warning('Skipping Subj %d, Stim %d: Load error.', subj, stim); dtl.clearData(); continue; end
        rawSignals = dtl.getRawSignals();
        currentMaxWindow = dtl.currentDuration; currentFs = dtl.fs; currentNChannels = dtl.nChannels;
        if ppc.fs ~= currentFs, ppc = PreProcessor('fs', currentFs); end
        if ordc.fs ~= currentFs, ordc = ORDCalculator(currentFs); end

        % 2. Preprocess Data
        processedSignals = rawSignals;
        if config.ppc_useZanoteli
            maxWinPPC = min(dtl.zanoteliSuggestedMMax(stim), size(processedSignals, 2));
            processedSignals = ppc.applyZanoteliPreprocessing(processedSignals, maxWinPPC);
        end
        if config.ppc_useFilter && ~isempty(processedSignals), processedSignals = ppc.applyAntunesFiltering(processedSignals); end
        if isempty(processedSignals) || size(processedSignals, 2) < 2, warning('Skipping Subj %d, Stim %d: No data after PPC.', subj, stim); dtl.clearData(); clear rawSignals processedSignals; continue; end
        currentMaxWindow = size(processedSignals, 2);

        % 3. Inner loop for ORD / Tester parameters
        for iStart = 1:numel(config.ord_startWindows)
             startWin = config.ord_startWindows(iStart);
             % --- Conditional Inner Loop ---
             for iParam = 1:numel(iterParam) % Iterate over K or stepSize
                 loopVal = iterParam(iParam); % Current K or stepSize value
             % ---
                 processedCount = processedCount + 1;

                 % --- ETA Calculation and Progress ---
                 if mod(processedCount, 50) == 0 || processedCount == totalIterations
                     elapsedTime = toc(startTime); progressStr = sprintf('  Processed %d / %d iterations...', processedCount, totalIterations);
                     if elapsedTime > 60 && processedCount > 0
                         timePerIter = elapsedTime / processedCount; remainingIter = totalIterations - processedCount; etaSeconds = remainingIter * timePerIter;
                         etaMinutes = floor(etaSeconds / 60); etaSecRem = round(mod(etaSeconds, 60)); etaStr = sprintf(' ETA: %d min %02d sec', etaMinutes, etaSecRem); progressStr = [progressStr, etaStr]; %#ok<AGROW>
                     end; fprintf('%s\n', progressStr);
                 end; % --- End ETA ---

                 if startWin >= currentMaxWindow, continue; end

                 % 4. Calculate MSC (Pass appropriate parameter)
                 mscOpts = {'method', config.ord_method, 'startWindow', startWin, 'maxWindow', currentMaxWindow};
                 if strcmp(config.ord_method, 'fixedK')
                     mscOpts = [mscOpts, {'K', loopVal}]; K_actual_for_key = loopVal; stepSize_for_key = NaN;
                 else % fixedStep
                     mscOpts = [mscOpts, {'stepSize', loopVal}]; stepSize_for_key = loopVal; K_actual_for_key = NaN;
                 end
                 [msc, epochs, M, K_actual] = ordc.calculateMSC(processedSignals, mscOpts{:});

                 if isempty(msc) || K_actual < 1 || K_actual>15 ||M < 6, continue; end

                 % --- Update Key and Params ---
                 if strcmp(config.ord_method, 'fixedK')
                     paramKey = sprintf('SW%d_K%d_Mthd%s', startWin, K_actual_for_key, config.ord_method);
                     currentOrdParams = struct('startWindow', startWin, 'K', K_actual_for_key, 'method', config.ord_method, 'M', M);
                 else % fixedStep
                     paramKey = sprintf('SW%d_SS%d_Mthd%s', startWin, stepSize_for_key, config.ord_method);
                     currentOrdParams = struct('startWindow', startWin, 'stepSize', stepSize_for_key, 'method', config.ord_method, 'M', M, 'K_actual', K_actual);
                 end
                 % ---

                 % 5. Run Tester
                 [decisionResults, ~, ~] = tester.runBetaCGST(msc, M, K_actual);

                 % 6. Extract counts and time
                 tp_matrix = decisionResults.TP(1:dtl.nSignalFrequencies, :, :); fn_matrix = decisionResults.FN(1:dtl.nSignalFrequencies, :, :);
                 fp_matrix = decisionResults.FP(dtl.nSignalFrequencies+1:end, :, :); tn_matrix = decisionResults.TN(dtl.nSignalFrequencies+1:end, :, :);
                 tp_this = sum(any(tp_matrix, 2), [1, 3]); fn_this = sum(any(fn_matrix, 2), [1, 3]);
                 fp_this = sum(any(fp_matrix, 2), [1, 3]); tn_this = sum(any(tn_matrix, 2), [1, 3]);
                 nSigOpps_this = dtl.nSignalFrequencies * currentNChannels; nNoiseOpps_this = dtl.nNoiseFrequencies * currentNChannels;

                 avg_det_time_this = NaN; n_det_events_this = 0;
                 if tp_this > 0 && ~isempty(epochs)
                     [~, stage_idx, ~] = ind2sub(size(tp_matrix), find(tp_matrix));
                     if ~isempty(stage_idx)
                         all_epoch_end_times = epochs(stage_idx, 2);
                         if ~isempty(all_epoch_end_times)
                             avg_det_time_this = mean(all_epoch_end_times); n_det_events_this = numel(all_epoch_end_times);
                         end
                     end
                 end

                 % 7. Aggregate into the map (SUMMING across subjects/stimuli)
                 if aggregatedResultsMap.isKey(paramKey)
                     aggData = aggregatedResultsMap(paramKey);
                     aggData.sumTP = aggData.sumTP + tp_this; aggData.sumFN = aggData.sumFN + fn_this;
                     aggData.sumFP = aggData.sumFP + fp_this; aggData.sumTN = aggData.sumTN + tn_this;
                     aggData.totalSignalOpps = aggData.totalSignalOpps + nSigOpps_this;
                     aggData.totalNoiseOpps = aggData.totalNoiseOpps + nNoiseOpps_this;
                     if ~isfield(aggData, 'MList'), aggData.MList = []; end
                     aggData.MList = [aggData.MList, M];
                     if ~isnan(avg_det_time_this) && n_det_events_this > 0
                         aggData.sumAvgDetTime = aggData.sumAvgDetTime + avg_det_time_this * n_det_events_this;
                         aggData.nAvgDetTime = aggData.nAvgDetTime + n_det_events_this;
                     end
                     aggregatedResultsMap(paramKey) = aggData;
                 else
                      aggregatedResultsMap(paramKey) = struct(...
                          'ordParams', currentOrdParams,... % Store ORD params + M
                          'sumTP', tp_this, 'sumFN', fn_this, 'sumFP', fp_this, 'sumTN', tn_this,...
                          'totalSignalOpps', nSigOpps_this, 'totalNoiseOpps', nNoiseOpps_this, ...
                          'MList', M, 'sumAvgDetTime', 0, 'nAvgDetTime', 0 ...
                          );
                      if ~isnan(avg_det_time_this) && n_det_events_this > 0
                           aggData = aggregatedResultsMap(paramKey);
                           aggData.sumAvgDetTime = avg_det_time_this * n_det_events_this;
                           aggData.nAvgDetTime = n_det_events_this;
                           aggregatedResultsMap(paramKey) = aggData;
                      end
                 end
                 clear msc decisionResults epochs tp_matrix fn_matrix fp_matrix tn_matrix;
             end % K or stepSize loop
        end % startWindow loop
        dtl.clearData();
        clear rawSignals processedSignals;
    end % Stimulus loop
end % Subject loop

if mod(processedCount, 50) ~= 0 && processedCount == totalIterations
     fprintf('  Processed %d / %d iterations...\n', processedCount, totalIterations);
end
fprintf('[%s] Experiment batch processing finished. Total time: %.2f seconds.\n', datetime, toc(startTime));

end