function aggregatedResultsMap = runSimulationBatch(config, dtl, ppc, ordc, tester)
% Runs the simulation loop for given configuration and objects.
% Returns map where key includes Sim Params + ORD Params, value has results for that SINGLE run.

fprintf('[%s] Starting simulation batch processing...\n', datetime);

aggregatedResultsMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
totalIterations = numel(config.sim_durations) * numel(config.sim_noiseMeans) * numel(config.sim_noiseStds) * numel(config.ord_startWindows) * numel(config.ord_K_stages);
processedCount = 0;

nChannels = dtl.nChannels;
nSignalFreqs = dtl.nSignalFrequencies;
nNoiseFreqs = dtl.nNoiseFrequencies;

startTime = tic; % Start timer for ETA calculation

for iDur = 1:numel(config.sim_durations)
    dur = config.sim_durations(iDur);
    for iMean = 1:numel(config.sim_noiseMeans)
        noiseMean = config.sim_noiseMeans(iMean);
        for iStd = 1:numel(config.sim_noiseStds)
            noiseStd = config.sim_noiseStds(iStd);

            % 1. Generate Simulated Data
            success = dtl.generateSimulation('duration', dur, 'noiseMean', noiseMean, 'noiseStd', noiseStd);
            if ~success, warning('Skipping Sim Dur=%d, Mean=%.1f, Std=%.1f: Gen error.', dur, noiseMean, noiseStd); dtl.clearData(); continue; end
            rawSignals = dtl.getRawSignals();
            currentMaxWindow = dtl.currentDuration;

            % 2. Preprocess Data
            processedSignals = rawSignals;
            if config.ppc_useZanoteli
                processedSignals = ppc.applyZanoteliPreprocessing(processedSignals, size(processedSignals, 2));
            end
            if config.ppc_useFilter && ~isempty(processedSignals)
                processedSignals = ppc.applyAntunesFiltering(processedSignals);
            end
             if isempty(processedSignals) || size(processedSignals, 2) < 2
                warning('Skipping Sim Dur=%d, Mean=%.1f, Std=%.1f: No data after PPC.', dur, noiseMean, noiseStd);
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
                     if mod(processedCount, 100) == 0 || processedCount == totalIterations
                         elapsedTime = toc(startTime);
                         progressStr = sprintf('  Processed %d / %d iterations...', processedCount, totalIterations);
                         if elapsedTime > 60 && processedCount > 0 % Show ETA after 60s
                             timePerIter = elapsedTime / processedCount;
                             remainingIter = totalIterations - processedCount;
                             etaSeconds = remainingIter * timePerIter;
                             % Format ETA
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

                     % --- Create Unique Key including Sim Params AND ORD Params ---
                     paramKey = sprintf('Dur%d_nM%.1f_nS%.1f_SW%d_K%d_Mthd%s', ...
                                        dur, noiseMean, noiseStd, startWin, K_val, config.ord_method);
                     currentFullParams = struct('startWindow', startWin, 'K', K_val, 'method', config.ord_method,...
                                                'duration', dur, 'noiseMean', noiseMean, 'noiseStd', noiseStd, 'M', M);

                     % 5. Run Tester
                     [decisionResults, ~, ~] = tester.runBetaCGST(msc, M, K_actual);

                     % 6. Extract counts and times
                     tp_matrix = decisionResults.TP(1:nSignalFreqs, :, :);
                     fn_matrix = decisionResults.FN(1:nSignalFreqs, :, :);
                     fp_matrix = decisionResults.FP(nSignalFreqs+1:end, :, :);
                     tn_matrix = decisionResults.TN(nSignalFreqs+1:end, :, :);
                     tp_this = sum(any(tp_matrix, 2), [1, 3]); fn_this = sum(any(fn_matrix, 2), [1, 3]);
                     fp_this = sum(any(fp_matrix, 2), [1, 3]); tn_this = sum(any(tn_matrix, 2), [1, 3]);
                     nSigOpps_this = nSignalFreqs * nChannels; nNoiseOpps_this = nNoiseFreqs * nChannels;

                     min_det_time_this = NaN; avg_det_time_this = NaN; n_det_events_this = 0;
                     if tp_this > 0 && ~isempty(epochs)
                         [~, stage_idx, ~] = ind2sub(size(tp_matrix), find(tp_matrix));
                         if ~isempty(stage_idx)
                             all_epoch_end_times = epochs(stage_idx, 2);
                             if ~isempty(all_epoch_end_times)
                                 min_det_time_this = min(all_epoch_end_times);
                                 avg_det_time_this = mean(all_epoch_end_times);
                                 n_det_events_this = numel(all_epoch_end_times);
                             end
                         end
                     end

                     % 7. Store results for THIS unique run (no aggregation here yet)
                      aggregatedResultsMap(paramKey) = struct(...
                          'fullParams', currentFullParams,...
                          'sumTP', tp_this, 'sumFN', fn_this, 'sumFP', fp_this, 'sumTN', tn_this,...
                          'totalSignalOpps', nSigOpps_this, 'totalNoiseOpps', nNoiseOpps_this, ...
                          'sumMinDetTime', min_det_time_this, 'nMinDetTime', double(~isnan(min_det_time_this)), ...
                          'sumAvgDetTime', avg_det_time_this * n_det_events_this, 'nAvgDetTime', n_det_events_this ...
                          );

                     clear msc decisionResults epochs tp_matrix fn_matrix fp_matrix tn_matrix;
                 end % K loop
            end % startWindow loop
            dtl.clearData();
            clear rawSignals processedSignals;
        end % noiseStd loop
    end % noiseMean loop
end % duration loop

% Final progress update if loop finished before hitting a multiple of 100
if mod(processedCount, 100) ~= 0 && processedCount == totalIterations
     fprintf('  Processed %d / %d iterations...\n', processedCount, totalIterations);
end

fprintf('[%s] Simulation batch processing finished. Total time: %.2f seconds.\n', datetime, toc(startTime));

end