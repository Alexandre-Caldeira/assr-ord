%% Clear workspace and initialize parameters
clearvars; close all; clc;
rng('default'); % For reproducible simulations

fprintf('[%s] Starting Refactored Simulation Script...\n', datetime);

%% --- Configuration ---
% Simulation Parameters (Loop Variables)
param_durations = [90];      % Simulation duration in seconds
param_noiseMeans = [-8 -8 0 -15]; % Mean SNR dB values to test
param_noiseStds = [1 0.1 1];       % Std Dev SNR dB values to test (can add more)

% Preprocessing
useAntunesFilter = false; % Set to true to apply filtering
useZanoteliPPC = true;    % Set to true to apply Zanoteli artifact rejection

% ORD Calculation Parameters (Loop Variables)
param_startWindows = 1:30; % Example: Reduced range for faster sim
param_K_stages = 2:5;   % Number of stages (tests) to try
param_method = 'fixedK'; % Use fixed number of stages ('fixedK' or 'fixedStep')

% ORDTester Parameters
tester_desiredAlpha = 0.05;

% Aggregation Variables
paramList = []; % Store parameters for each result entry
aggregated_TP = [];
aggregated_TN = [];
aggregated_FP = [];
aggregated_FN = [];
aggregated_MinDetectTime = []; % Store minimum detection time (TP only)
aggregated_AvgDetectTime = []; % Store average detection time (TP only)

nTotalSignalFreqTests = 0;
nTotalNoiseFreqTests = 0;

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
nChannels = dtl.nChannels;

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
            simParams = struct('duration', dur, 'noiseMean', noiseMean, 'noiseStd', noiseStd);
            success = dtl.generateSimulation('duration', dur, 'noiseMean', noiseMean, 'noiseStd', noiseStd);
            if ~success
                warning('Skipping Sim Dur=%d, Mean=%.1f, Std=%.1f due to generation error.', dur, noiseMean, noiseStd);
                dtl.clearData();
                continue;
            end
            rawSignals = dtl.getRawSignals();
            currentMaxWindow = dtl.currentDuration;

            % 2. Preprocess Data (Optional for Sim, but good practice)
             if useZanoteliPPC
                % For sim, no suggested max M, use all available
                maxWindowsForPPC = size(rawSignals, 2); 
                processedSignals = ppc.applyZanoteliPreprocessing(rawSignals, maxWindowsForPPC);
             else
                 processedSignals = rawSignals; % Skip Zanoteli
             end
             
             if useAntunesFilter
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

                     tp_count = sum(any(tp_matrix, [2, 3]));
                     fn_count = sum(any(fn_matrix, [2, 3]));
                     fp_count = sum(any(fp_matrix, [2, 3]));
                     tn_count = sum(any(tn_matrix, [2, 3]));

                     % Calculate Detection Times (only for TPs)
                     min_det_time = NaN;
                     avg_det_time = NaN;
                     if tp_count > 0
                         detection_times = [];
                         % Find the stage index (k) where TP first occurs for each freq/channel
                         [freq_idx, stage_idx, chan_idx] = ind2sub(size(tp_matrix), find(tp_matrix));
                         
                         % Get the *end* window index for the stage where detection occurred
                         % epochs is K x 2, where epochs(k, 2) is the end window index for stage k
                         if ~isempty(stage_idx)
                             epoch_end_times = epochs(stage_idx, 2); % Time in seconds (assuming 1s windows)
                             detection_times = epoch_end_times;
                             
                             min_det_time = min(detection_times);
                             avg_det_time = mean(detection_times);
                         end
                     end

                     % Store aggregated counts and times for this parameter combination
                     paramList = [paramList; currentParams]; %#ok<AGROW>
                     aggregated_TP = [aggregated_TP; tp_count]; %#ok<AGROW>
                     aggregated_FN = [aggregated_FN; fn_count]; %#ok<AGROW>
                     aggregated_FP = [aggregated_FP; fp_count]; %#ok<AGROW>
                     aggregated_TN = [aggregated_TN; tn_count]; %#ok<AGROW>
                     aggregated_MinDetectTime = [aggregated_MinDetectTime; min_det_time]; %#ok<AGROW>
                     aggregated_AvgDetectTime = [aggregated_AvgDetectTime; avg_det_time]; %#ok<AGROW>

                     nTotalSignalFreqTests = nTotalSignalFreqTests + nSignalFreqs;
                     nTotalNoiseFreqTests = nTotalNoiseFreqTests + nNoiseFreqs;

                     clear msc decisionResults epochs tp_matrix fn_matrix fp_matrix tn_matrix;

                 end % K loop
            end % startWindow loop

            dtl.clearData();
            clear rawSignals processedSignals;
            fprintf('[%s] Finished Sim Dur=%ds, SNR Mean=%.1f, SNR Std=%.1f.\n', datetime, dur, noiseMean, noiseStd);

        end % noiseStd loop
    end % noiseMean loop
end % duration loop

fprintf('\n[%s] All simulation configurations processed.\n', datetime);

%% --- Post-Processing and Visualization ---

if isempty(paramList)
    error('No results were generated. Check simulation parameters.');
end

% Calculate Rates
tp_rate = aggregated_TP / nSignalFreqs;
fn_rate = aggregated_FN / nSignalFreqs;
fp_rate = aggregated_FP / nNoiseFreqs; % Handle nNoiseFreqs = 0
tn_rate = aggregated_TN / nNoiseFreqs; % Handle nNoiseFreqs = 0
if nNoiseFreqs == 0
    fp_rate(:) = NaN;
    tn_rate(:) = NaN;
end

% Create results table
resultsTable = struct2table(paramList);
resultsTable.TP_Rate = tp_rate * 100;
resultsTable.FN_Rate = fn_rate * 100;
resultsTable.FP_Rate = fp_rate * 100;
resultsTable.TN_Rate = tn_rate * 100;
resultsTable.MinDetectTime_s = aggregated_MinDetectTime;
resultsTable.AvgDetectTime_s = aggregated_AvgDetectTime;

disp('Simulation Results Summary Table (Rates in %):');
disp(head(resultsTable, 10)); % Show first 10 rows

% Overall Average Rates
avg_tp_rate = mean(resultsTable.TP_Rate, 'omitnan');
avg_fn_rate = mean(resultsTable.FN_Rate, 'omitnan');
avg_fp_rate = mean(resultsTable.FP_Rate, 'omitnan');
avg_tn_rate = mean(resultsTable.TN_Rate, 'omitnan');
avg_min_time = mean(resultsTable.MinDetectTime_s, 'omitnan');
avg_avg_time = mean(resultsTable.AvgDetectTime_s, 'omitnan');


confmat_avg = table(avg_fn_rate, avg_fp_rate, avg_tp_rate, avg_tn_rate, avg_min_time, avg_avg_time, ...
    'VariableNames', {'Avg FN Rate (%)', 'Avg FP Rate (%)', 'Avg TP Rate (%)', 'Avg TN Rate (%)', 'Avg Min Detect Time (s)', 'Avg Avg Detect Time (s)'});
disp('Overall Average Confusion Rates (%) and Detection Times (s):');
disp(confmat_avg);

% Example Plots (Adapt as needed)
figure;
subplot(2, 2, 1); boxchart(resultsTable.TP_Rate); title('TP Rate (%)'); grid on; ylim([0 100]);
subplot(2, 2, 2); boxchart(resultsTable.FP_Rate); title('FP Rate (%)'); grid on; ylim([0 100]);
subplot(2, 2, 3); boxchart(resultsTable.FN_Rate); title('FN Rate (%)'); grid on; ylim([0 100]);
subplot(2, 2, 4); boxchart(resultsTable.TN_Rate); title('TN Rate (%)'); grid on; ylim([0 100]);
sgtitle('Distribution of Rates Across All Simulation Parameter Sets');

figure;
subplot(1, 2, 1); boxchart(resultsTable.MinDetectTime_s); title('Min Detection Time (s)'); grid on;
subplot(1, 2, 2); boxchart(resultsTable.AvgDetectTime_s); title('Avg Detection Time (s)'); grid on;
sgtitle('Distribution of Detection Times (TP only) Across All Simulation Parameter Sets');


% Pareto Front Plot (TP Rate vs Min Detection Time)
figure;
valid_idx = ~isnan(resultsTable.MinDetectTime_s) & resultsTable.TP_Rate > 0;
if sum(valid_idx) > 1
    tp_vals = resultsTable.TP_Rate(valid_idx);
    time_vals = resultsTable.MinDetectTime_s(valid_idx);

    % Pareto front expects minimization, so minimize (-TP_Rate) and Time
    [p, p_indices] = paretoFront([-tp_vals, time_vals]);

    % Sort pareto points by time for plotting
    [sorted_time, sort_idx] = sort(p(:, 2));
    sorted_neg_tp = p(sort_idx, 1);

    plot(time_vals, tp_vals, 'k.', 'MarkerSize', 8, 'DisplayName', 'All Configurations');
    hold on;
    plot(sorted_time, -sorted_neg_tp, '-ob', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Pareto Front');
    xlabel('Minimum Detection Time (s)');
    ylabel('True Positive Rate (%)');
    title('Pareto Front: TP Rate vs. Minimum Detection Time');
    legend('Location', 'SouthEast');
    grid on;
    ylim([0 105]); % Ensure full range is visible
else
    warning('Not enough valid data points to generate Pareto plot.');
end


% --- Save Results ---
resultsFilename = sprintf('sim_refactored_results_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
try
    save(resultsFilename, 'resultsTable', 'confmat_avg', 'param_durations', 'param_noiseMeans', 'param_noiseStds', 'param_startWindows', 'param_K_stages');
    fprintf('[%s] Simulation results saved to %s\n', datetime, resultsFilename);
catch ME_save
    warning('Could not save results: %s', ME_save.message);
end

fprintf('[%s] Simulation script finished.\n', datetime);