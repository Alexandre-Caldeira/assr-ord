function aggregatedResultsMap = runSimulationBatch(config, dtl, ppc, ordc, tester)
% Runs the simulation loop for given configuration and objects.
% Returns a map where keys are ORD param strings and values are structs
% containing SUMMED counts/times/opportunities across simulation runs.

fprintf('[%s] Starting simulation batch processing...\n', datetime);

aggregatedResultsMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
paramIdx = 0;
totalIterations = numel(config.sim_durations) * numel(config.sim_noiseMeans) * numel(config.sim_noiseStds) * numel(config.ord_startWindows) * numel(config.ord_K_stages);
processedCount = 0;

% Get fixed N channels for simulation
nChannels = dtl.nChannels;
nSignalFreqs = dtl.nSignalFrequencies;
nNoiseFreqs = dtl.nNoiseFrequencies;

for iDur = 1:numel(config.sim_durations)
    dur = config.sim_durations(iDur);
    for iMean = 1:numel(config.sim_noiseMeans)
        noiseMean = config.sim_noiseMeans(iMean);
        for iStd = 1:numel(config.sim_noiseStds)
            noiseStd = config.sim_noiseStds(iStd);

            % fprintf('\n[%s] Sim: Dur=%ds, SNR Mean=%.1f, Std=%.1f...\n', datetime, dur, noiseMean, noiseStd); % Verbose

            % 1. Generate Simulated Data
            success = dtl.generateSimulation('duration', dur, 'noiseMean', noiseMean, 'noiseStd', noiseStd);
            if ~success, warning('Skipping Sim Dur=%d, Mean=%.1f, Std=%.1f: Gen error.', dur, noiseMean, noiseStd); dtl.clearData(); continue; end
            rawSignals = dtl.getRawSignals();
            currentMaxWindow = dtl.currentDuration;

            % 2. Preprocess Data (Optional for Sim)
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
                     paramIdx = paramIdx + 1;
                     processedCount = processedCount + 1;
                     if mod(processedCount, 100) == 0 % Progress marker
                         fprintf('  Processed %d / %d iterations...\n', processedCount, totalIterations);
                     end

                     paramKey = sprintf('SW%d_K%d_Mthd%s', startWin, K_val, config.ord_method);
                     currentOrdParams = struct('startWindow', startWin, 'K', K_val, 'method', config.ord_method);

                     if startWin >= currentMaxWindow, continue; end

                     % 4. Calculate MSC
                     [msc, epochs, M, K_actual] = ordc.calculateMSC(processedSignals, ...
                                                                     'method', config.ord_method, ...
                                                                     'startWindow', startWin, ...
                                                                     'K', K_val, ...
                                                                     'maxWindow', currentMaxWindow);

                     if isempty(msc) || K_actual < 1 || M < 2, continue; end

                     % 5. Run Tester
                     [decisionResults, ~, ~] = tester.runBetaCGST(msc, M, K_actual);

                     % 6. Extract counts and times for THIS sim run / ORD param config
                     tp_matrix = decisionResults.TP(1:nSignalFreqs, :, :);
                     fn_matrix = decisionResults.FN(1:nSignalFreqs, :, :);
                     fp_matrix = decisionResults.FP(nSignalFreqs+1:end, :, :);
                     tn_matrix = decisionResults.TN(nSignalFreqs+1:end, :, :);

                     tp_this = sum(any(tp_matrix, 2), [1, 3]);
                     fn_this = sum(any(fn_matrix, 2), [1, 3]);
                     fp_this = sum(any(fp_matrix, 2), [1, 3]);
                     tn_this = sum(any(tn_matrix, 2), [1, 3]);

                     nSigOpps_this = nSignalFreqs * nChannels;
                     nNoiseOpps_this = nNoiseFreqs * nChannels;

                     % Calculate Detection Times (only for TPs)
                     min_det_time_this = NaN;
                     avg_det_time_this = NaN;
                     n_det_events_this = 0;
                     if tp_this > 0 && ~isempty(epochs)
                         [~, stage_idx, ~] = ind2sub(size(tp_matrix), find(tp_matrix));
                         if ~isempty(stage_idx)
                             all_epoch_end_times = epochs(stage_idx, 2); % Time = end window index
                             if ~isempty(all_epoch_end_times)
                                 min_det_time_this = min(all_epoch_end_times);
                                 avg_det_time_this = mean(all_epoch_end_times); % Mean time for this config
                                 n_det_events_this = numel(all_epoch_end_times); % How many TPs contributed
                             end
                         end
                     end

                     % 7. Aggregate into the map
                     if aggregatedResultsMap.isKey(paramKey)
                         aggData = aggregatedResultsMap(paramKey);
                         aggData.sumTP = aggData.sumTP + tp_this;
                         aggData.sumFN = aggData.sumFN + fn_this;
                         aggData.sumFP = aggData.sumFP + fp_this;
                         aggData.sumTN = aggData.sumTN + tn_this;
                         aggData.totalSignalOpps = aggData.totalSignalOpps + nSigOpps_this;
                         aggData.totalNoiseOpps = aggData.totalNoiseOpps + nNoiseOpps_this;
                         % Aggregate times carefully (handle NaNs)
                         if ~isnan(min_det_time_this)
                             aggData.sumMinDetTime = aggData.sumMinDetTime + min_det_time_this;
                             aggData.nMinDetTime = aggData.nMinDetTime + 1;
                         end
                          if ~isnan(avg_det_time_this) && n_det_events_this > 0
                             aggData.sumAvgDetTime = aggData.sumAvgDetTime + avg_det_time_this * n_det_events_this; % Sum of times
                             aggData.nAvgDetTime = aggData.nAvgDetTime + n_det_events_this; % Total detections
                         end
                         aggregatedResultsMap(paramKey) = aggData;
                     else
                         % Create new entry, initialize time fields
                         aggregatedResultsMap(paramKey) = struct(...
                              'ordParams', currentOrdParams,...
                              'sumTP', tp_this, 'sumFN', fn_this, 'sumFP', fp_this, 'sumTN', tn_this,...
                              'totalSignalOpps', nSigOpps_this, 'totalNoiseOpps', nNoiseOpps_this, ...
                              'sumMinDetTime', 0, 'nMinDetTime', 0, ... % Init time accumulators
                              'sumAvgDetTime', 0, 'nAvgDetTime', 0 ...
                              );
                          % Add first time value if valid
                          if ~isnan(min_det_time_this)
                             aggData = aggregatedResultsMap(paramKey);
                             aggData.sumMinDetTime = min_det_time_this;
                             aggData.nMinDetTime = 1;
                             aggregatedResultsMap(paramKey) = aggData;
                          end
                           if ~isnan(avg_det_time_this) && n_det_events_this > 0
                             aggData = aggregatedResultsMap(paramKey);
                             aggData.sumAvgDetTime = avg_det_time_this * n_det_events_this;
                             aggData.nAvgDetTime = n_det_events_this;
                             aggregatedResultsMap(paramKey) = aggData;
                          end
                     end

                     clear msc decisionResults epochs tp_matrix fn_matrix fp_matrix tn_matrix;
                 end % K loop
            end % startWindow loop

            dtl.clearData();
            clear rawSignals processedSignals;
            % fprintf('[%s] Finished Sim Dur=%ds, SNR Mean=%.1f, Std=%.1f.\n', datetime, dur, noiseMean, noiseStd); % Verbose
        end % noiseStd loop
    end % noiseMean loop
end % duration loop

fprintf('[%s] Simulation batch processing finished. Processed %d iterations.\n', datetime, processedCount);

end % End function runSimulationBatch