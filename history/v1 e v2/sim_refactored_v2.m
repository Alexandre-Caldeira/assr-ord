%% Clear workspace and initialize parameters
clearvars; close all; clc;
rng('default'); % For reproducible simulations

fprintf('[%s] Starting Refactored Simulation Script...\n', datetime);

%% --- Configuration ---
% Simulation Parameters (Loop Variables)
param_durations = [30];      % Simulation duration in seconds
param_noiseMeans = [-30 -20]; % Mean SNR dB values to test
param_noiseStds = [5];       % Std Dev SNR dB values to test (can add more)

% Preprocessing
useAntunesFilter = false; % Set to true to apply filtering
useZanoteliPPC = true;    % Set to true to apply Zanoteli artifact rejection (less common for sim, but possible)

% ORD Calculation Parameters (Loop Variables)
param_startWindows = 1:5; % Example: Reduced range for faster sim
param_K_stages = 2:5;   % Number of stages (tests) to try
param_method = 'fixedK'; % Use fixed number of stages ('fixedK' or 'fixedStep')

% ORDTester Parameters
tester_desiredAlpha = 0.05;

% Aggregation Variables
paramList = []; % Store parameters for each result entry
aggregated_TP_counts = [];
aggregated_FN_counts = [];
aggregated_FP_counts = [];
aggregated_TN_counts = [];
aggregated_nSignalOpps = []; % Store opportunities per config
aggregated_nNoiseOpps = [];
aggregated_MinDetectTime = []; % Store minimum detection time (TP only)
aggregated_AvgDetectTime = []; % Store average detection time (TP only)


%% --- Initialization ---
try
    % Initialize DataLoader in 'sim' mode with placeholder values
    dtl = DataLoader('sim', 'simDuration', param_durations(1), ...
                     'noiseMean', param_noiseMeans(1), 'noiseStd', param_noiseStds(1));
    ppc = PreProcessor('fs', dtl.fs);
    ordc = ORDCalculator(dtl.fs);
    tester = ORDTester(dtl.signalFrequencies, dtl.noiseFrequencies, 'desiredAlpha', tester_desiredAlpha);
catch ME
    error('Initialization failed: %s', ME.message);
end

nSignalFreqs = tester.nSignalFrequencies;
nNoiseFreqs = tester.nNoiseFrequencies;
nChannels = dtl.nChannels; % Fs and nChannels are fixed in sim mode usually

%% --- Parameter Iteration Loop ---
paramIdx = 0;
totalConfigurations = numel(param_durations) * numel(param_noiseMeans) * numel(param_noiseStds) * numel(param_startWindows) * numel(param_K_stages);
fprintf('[%s] Total simulation configurations to process: %d\n', datetime, totalConfigurations);

for iDur = 1:numel(param_durations)
    dur = param_durations(iDur);
    for iMean = 1:numel(param_noiseMeans)
        noiseMean = param_noiseMeans(iMean);
        for iStd = 1:numel(param_noiseStds)
            noiseStd = param_noiseStds(iStd);

            fprintf('\n[%s] Processing Sim: Dur=%ds, SNR Mean=%.1f, SNR Std=%.1f...\n', datetime, dur, noiseMean, noiseStd);

            % 1. Generate Simulated Data
            success = dtl.generateSimulation('duration', dur, 'noiseMean', noiseMean, 'noiseStd', noiseStd);
            if ~success
                warning('Skipping Sim Dur=%d, Mean=%.1f, Std=%.1f due to generation error.', dur, noiseMean, noiseStd);
                dtl.clearData();
                continue;
            end
            rawSignals = dtl.getRawSignals();
            currentMaxWindow = dtl.currentDuration;

            % 2. Preprocess Data (Optional for Sim)
            processedSignals = rawSignals; % Start with raw
            if useZanoteliPPC
                maxWindowsForPPC = size(processedSignals, 2); % Use all available for sim
                processedSignals = ppc.applyZanoteliPreprocessing(processedSignals, maxWindowsForPPC);
            end
            if useAntunesFilter && ~isempty(processedSignals)
                processedSignals = ppc.applyAntunesFiltering(processedSignals);
            end

             if isempty(processedSignals) || size(processedSignals, 2) < 2
                warning('Skipping Sim Dur=%d, Mean=%.1f, Std=%.1f: Not enough data after preprocessing (%d windows).', ...
                        dur, noiseMean, noiseStd, size(processedSignals, 2));
                dtl.clearData();
                clear rawSignals processedSignals;
                continue;
             end
             currentMaxWindow = size(processedSignals, 2); % Update max window

            % --- Inner loop for ORD / Tester parameters ---
            for iStart = 1:numel(param_startWindows)
                 startWin = param_startWindows(iStart);
                 for iK = 1:numel(param_K_stages)
                     K_val = param_K_stages(iK);
                     paramIdx = paramIdx + 1;
                     fprintf('--- Param Set %d/%d: startWin=%d, K=%d ---\n', paramIdx, totalConfigurations, startWin, K_val);

                     currentParams = struct('duration', dur, 'noiseMean', noiseMean, 'noiseStd', noiseStd, ...
                                            'startWindow', startWin, 'K', K_val, 'method', param_method);

                     if startWin >= currentMaxWindow
                          warning('Skipping param set: startWindow (%d) >= maxWindow (%d) for this sim.', startWin, currentMaxWindow);
                          continue;
                     end

                     % 3. Calculate MSC
                     [msc, epochs, M, K_actual] = ordc.calculateMSC(processedSignals, ...
                                                                     'method', param_method, ...
                                                                     'startWindow', startWin, ...
                                                                     'K', K_val, ...
                                                                     'maxWindow', currentMaxWindow);

                     if isempty(msc) || K_actual < 1 || M < 2
                         warning('Skipping param set: MSC calculation failed or resulted in invalid M/K (M=%d, K=%d).', M, K_actual);
                         continue;
                     end

                     % 4. Run Tester
                     [decisionResults, ~, ~] = tester.runBetaCGST(msc, M, K_actual);

                     % 5. Aggregate Results
                     tp_matrix = decisionResults.TP(1:nSignalFreqs, :, :); % SigFreqs x K x Channels
                     fn_matrix = decisionResults.FN(1:nSignalFreqs, :, :);
                     fp_matrix = decisionResults.FP(nSignalFreqs+1:end, :, :); % NoiseFreqs x K x Channels
                     tn_matrix = decisionResults.TN(nSignalFreqs+1:end, :, :);

                     tp_this_config = sum(any(tp_matrix, 2), [1, 3]);
                     fn_this_config = sum(any(fn_matrix, 2), [1, 3]);
                     fp_this_config = sum(any(fp_matrix, 2), [1, 3]);
                     tn_this_config = sum(any(tn_matrix, 2), [1, 3]);

                     % Calculate Detection Times (only for TPs)
                     min_det_time = NaN;
                     avg_det_time = NaN;
                     if tp_this_config > 0 && ~isempty(epochs) % Ensure epochs were calculated
                         detection_times = [];
                         % Find the stage index (k) where TP first occurs for each freq/channel pair
                         [~, stage_idx, ~] = ind2sub(size(tp_matrix), find(tp_matrix)); % Get stage index k

                         if ~isempty(stage_idx)
                             % Get the *end* window index for the stage where detection occurred
                             % epochs is K x 2, where epochs(k, 2) is the end window index for stage k
                             % Assuming 1s windows, end window index = time in seconds
                             unique_stage_indices = unique(stage_idx); % Get unique stages with detections
                             epoch_end_times = epochs(unique_stage_indices, 2); % Time in seconds
                             % We need to consider *all* detection events for average,
                             % but only the earliest overall for minimum.
                             all_epoch_end_times = epochs(stage_idx, 2); % Get end time for every TP event

                             if ~isempty(all_epoch_end_times)
                                 min_det_time = min(all_epoch_end_times);
                                 avg_det_time = mean(all_epoch_end_times);
                             end
                         end
                     end

                     % Store aggregated counts and times for this parameter combination
                     paramList = [paramList; currentParams]; %#ok<AGROW>
                     aggregated_TP_counts = [aggregated_TP_counts; tp_this_config]; %#ok<AGROW>
                     aggregated_FN_counts = [aggregated_FN_counts; fn_this_config]; %#ok<AGROW>
                     aggregated_FP_counts = [aggregated_FP_counts; fp_this_config]; %#ok<AGROW>
                     aggregated_TN_counts = [aggregated_TN_counts; tn_this_config]; %#ok<AGROW>
                     aggregated_MinDetectTime = [aggregated_MinDetectTime; min_det_time]; %#ok<AGROW>
                     aggregated_AvgDetectTime = [aggregated_AvgDetectTime; avg_det_time]; %#ok<AGROW>


                     % Denominators for rates for THIS configuration
                     nSignalOpportunities = nSignalFreqs * nChannels;
                     nNoiseOpportunities = nNoiseFreqs * nChannels;
                     aggregated_nSignalOpps = [aggregated_nSignalOpps; nSignalOpportunities]; %#ok<AGROW>
                     aggregated_nNoiseOpps = [aggregated_nNoiseOpps; nNoiseOpportunities]; %#ok<AGROW>


                     clear msc decisionResults epochs tp_matrix fn_matrix fp_matrix tn_matrix;

                 end % K loop
            end % startWindow loop

            dtl.clearData();
            clear rawSignals processedSignals;
            fprintf('[%s] Finished Sim Dur=%ds, SNR Mean=%.1f, SNR Std=%.1f.\n', datetime, dur, noiseMean, noiseStd);

        end % noiseStd loop
    end % noiseMean loop
end % duration loop

fprintf('\n[%s] All simulation configurations processed. Calculating final metrics...\n', datetime);

%% --- Post-Processing and Visualization ---

if isempty(paramList)
    error('No results were generated. Check simulation parameters.');
end

% Prepare aggregated data structure for MetricsCalculator
aggregatedData = struct(...
    'TP', aggregated_TP_counts, ...
    'FN', aggregated_FN_counts, ...
    'FP', aggregated_FP_counts, ...
    'TN', aggregated_TN_counts, ...
    'nSignalOpps', aggregated_nSignalOpps, ...
    'nNoiseOpps', aggregated_nNoiseOpps, ...
    'MinDetectTime', aggregated_MinDetectTime, ... % Add time data
    'AvgDetectTime', aggregated_AvgDetectTime ...
);

% Create MetricsCalculator instance
metrics = MetricsCalculator(paramList, aggregatedData, nSignalFreqs, nNoiseFreqs, nChannels);

% Display Summary
metrics.displaySummary();

% Generate Plots
metrics.plotRateDistributions(figure(1));
metrics.plotTimeDistributions(figure(2));
metrics.plotParetoFront('MinDetectTime_s', 'TP_Rate', figure(3)); % Pareto Time vs TP
metrics.plotParetoFront('FP_Rate', 'TP_Rate', figure(4)); % Pareto FP vs TP


% --- Save Results ---
resultsTable = metrics.getResultsTable(); % Get table for saving
confmat_avg = metrics.getAverageResults(); % Get averages for saving
resultsFilename = sprintf('sim_refactored_results_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
try
    % Save relevant config parameters along with results
    save(resultsFilename, 'resultsTable', 'confmat_avg', ...
         'param_durations', 'param_noiseMeans', 'param_noiseStds', ...
         'param_startWindows', 'param_K_stages', ...
         'nSignalFreqs', 'nNoiseFreqs', 'nChannels', 'tester_desiredAlpha', ...
         'useAntunesFilter', 'useZanoteliPPC');
    fprintf('[%s] Simulation results saved to %s\n', datetime, resultsFilename);
catch ME_save
    warning('Could not save results: %s', ME_save.message);
end

fprintf('[%s] Simulation script finished.\n', datetime);