%% Clear workspace and initialize parameters
clearvars; close all; clc;
rng('default'); % For reproducible noise frequency generation if needed

fprintf('[%s] Starting Refactored Experiment Script...\n', datetime);

%% --- Configuration ---
% Data Loading
expDataPath = 'C:\PPGEE\SBEB_CBA_24\CGST_figuras\Sinais_EEG\' %'C:\Users\alexa\Desktop\Sinais_EEG\'; % ADJUST AS NEEDED
selectedSubjects = 1:11;
selectedStimuli = 3; % Example: 50dB stimulus
useAntunesFilter = false; % Set to true to apply filtering

% ORD Calculation Parameters (Loop Variables)
% Example: Iterate through different K values for a fixed start window
param_startWindows = 5:140; % Fixed start window
param_K_stages = 2:10;   % Number of stages (tests) to try
param_method = 'fixedK'; % Use fixed number of stages ('fixedK' or 'fixedStep')
% param_stepSize = 24; % Only used if method = 'fixedStep'

% ORDTester Parameters
tester_desiredAlpha = 0.05;

% Aggregation Variables
results = struct('params', [], 'totalTP', 0, 'totalTN', 0, 'totalFP', 0, 'totalFN', 0, 'count', 0);
paramList = []; % Store parameters for each result entry
aggregated_TP = []; % Store TP count per parameter set
aggregated_TN = [];
aggregated_FP = [];
aggregated_FN = [];
nTotalSignalFreqTests = 0; % Denominator for rates
nTotalNoiseFreqTests = 0;

%% --- Initialization ---
try
    dtl = DataLoader('exp', 'dataPath', expDataPath);
    ppc = PreProcessor('fs', dtl.fs); % Match PreProcessor Fs to DataLoader Fs
    ordc = ORDCalculator(dtl.fs);
    tester = ORDTester(dtl.signalFrequencies, dtl.noiseFrequencies, 'desiredAlpha', tester_desiredAlpha);
catch ME
    error('Initialization failed: %s', ME.message);
end

nSignalFreqs = tester.nSignalFrequencies;
nNoiseFreqs = tester.nNoiseFrequencies;
nChannels = dtl.nChannels; % Assuming using all default channels

%% --- Parameter Iteration Loop ---
paramIdx = 0;
totalConfigurations = numel(selectedSubjects) * numel(selectedStimuli) * numel(param_startWindows) * numel(param_K_stages);
fprintf('[%s] Total configurations to process: %d\n', datetime, totalConfigurations);

for iSub = 1:numel(selectedSubjects)
    subj = selectedSubjects(iSub);
    for iStim = 1:numel(selectedStimuli)
        stim = selectedStimuli(iStim);

        fprintf('\n[%s] Processing Subject %d, Stimulus %d...\n', datetime, subj, stim);

        % 1. Load Data
        success = dtl.loadExperiment(subj, stim);
        if ~success
            warning('Skipping Subject %d, Stimulus %d due to loading error.', subj, stim);
            dtl.clearData(); % Clear any partial data
            continue; % Skip to next stimulus/subject
        end
        rawSignals = dtl.getRawSignals();
        currentMaxWindow = dtl.currentDuration; % Max available windows for this exam

        % 2. Preprocess Data
        % Use Zanoteli suggested max M, or just use all available after rejection
        maxWindowsForPPC = dtl.zanoteliSuggestedMMax(stim);
        % Ensure maxWindowsForPPC doesn't exceed actual available windows
        maxWindowsForPPC = min(maxWindowsForPPC, size(rawSignals, 2)); 
        
        processedSignals = ppc.applyZanoteliPreprocessing(rawSignals, maxWindowsForPPC);
        if useAntunesFilter
            processedSignals = ppc.applyAntunesFiltering(processedSignals);
        end
        
        % Check if preprocessing left enough data
        if isempty(processedSignals) || size(processedSignals, 2) < 2 
            warning('Skipping Subject %d, Stimulus %d: Not enough data after preprocessing (%d windows).', ...
                    subj, stim, size(processedSignals, 2));
            dtl.clearData();
            clear rawSignals processedSignals;
            continue;
        end
        currentMaxWindow = size(processedSignals, 2); % Update max window after preprocessing


        % --- Inner loop for ORD / Tester parameters ---
        for iStart = 1:numel(param_startWindows)
             startWin = param_startWindows(iStart);
             for iK = 1:numel(param_K_stages)
                 K_val = param_K_stages(iK);
                 paramIdx = paramIdx + 1;
                 fprintf('--- Param Set %d/%d: startWin=%d, K=%d ---\n', paramIdx, totalConfigurations, startWin, K_val);

                 currentParams = struct('subject', subj, 'stimulus', stim, ...
                                        'startWindow', startWin, 'K', K_val, 'method', param_method); % Add other params if they vary

                 % Check if start window is valid for current data
                 if startWin >= currentMaxWindow
                      warning('Skipping param set: startWindow (%d) >= maxWindow (%d) for this exam.', startWin, currentMaxWindow);
                      continue;
                 end

                 % 3. Calculate MSC
                 % ordParams = struct('method', param_method, 'startWindow', startWin, ...
                 %                    'K', K_val, 'maxWindow', currentMaxWindow);
                 % If using fixedStep, add stepSize: ordParams.stepSize = param_stepSize;
                 
                 [msc, epochs, M, K_actual] = ordc.calculateMSC(processedSignals, ...
                                                                 'method', param_method, ...
                                                                 'startWindow', startWin, ...
                                                                 'K', K_val, ...
                                                                 'maxWindow', currentMaxWindow);

                 if isempty(msc) || K_actual < 1 || M < 2
                     warning('Skipping param set: MSC calculation failed or resulted in invalid M/K (M=%d, K=%d).', M, K_actual);
                     continue; % Skip to next parameter set
                 end

                 % 4. Run Tester
                 [decisionResults, ~, ~] = tester.runBetaCGST(msc, M, K_actual);

                 % 5. Aggregate Results for this parameter set across channels
                 % Sum outcomes where decision was made (TP=1 or FN=1 or FP=1 or TN=1)
                 % Sum over stages (dim 2) and channels (dim 3)
                 
                 % TP: Sum over signal freqs where TP=1 at any stage/channel
                 tp_count = sum(any(decisionResults.TP(1:nSignalFreqs, :, :), [2, 3])); 
                 % FN: Sum over signal freqs where FN=1 at any stage/channel
                 fn_count = sum(any(decisionResults.FN(1:nSignalFreqs, :, :), [2, 3])); 
                 % FP: Sum over noise freqs where FP=1 at any stage/channel
                 fp_count = sum(any(decisionResults.FP(nSignalFreqs+1:end, :, :), [2, 3]));
                 % TN: Sum over noise freqs where TN=1 at any stage/channel
                 tn_count = sum(any(decisionResults.TN(nSignalFreqs+1:end, :, :), [2, 3]));
                 
                 % Store aggregated counts for this specific parameter combination
                 paramList = [paramList; currentParams]; %#ok<AGROW>
                 aggregated_TP = [aggregated_TP; tp_count]; %#ok<AGROW>
                 aggregated_FN = [aggregated_FN; fn_count]; %#ok<AGROW>
                 aggregated_FP = [aggregated_FP; fp_count]; %#ok<AGROW>
                 aggregated_TN = [aggregated_TN; tn_count]; %#ok<AGROW>
                 
                 % Update overall denominators (count each potential test once per config)
                 nTotalSignalFreqTests = nTotalSignalFreqTests + nSignalFreqs;
                 nTotalNoiseFreqTests = nTotalNoiseFreqTests + nNoiseFreqs;


                 % Clear large intermediate data for this param set if needed
                 clear msc decisionResults epochs;

             end % K loop
        end % startWindow loop

        % Clear data for the current subject/stimulus before loading next
        dtl.clearData();
        clear rawSignals processedSignals;
        fprintf('[%s] Finished Subject %d, Stimulus %d.\n', datetime, subj, stim);

    end % Stimulus loop
end % Subject loop

fprintf('\n[%s] All configurations processed.\n', datetime);

%% --- Post-Processing and Visualization ---

if isempty(paramList)
    error('No results were generated. Check parameters and data.');
end

% Calculate Rates for each parameter set
% Denominator should reflect the number of opportunities
% For TP/FN rate: number of signal frequencies tested in that config (nSignalFreqs)
% For FP/TN rate: number of noise frequencies tested in that config (nNoiseFreqs)
tp_rate = aggregated_TP / nSignalFreqs;
fn_rate = aggregated_FN / nSignalFreqs;
fp_rate = aggregated_FP / nNoiseFreqs; % Handle case nNoiseFreqs = 0
tn_rate = aggregated_TN / nNoiseFreqs; % Handle case nNoiseFreqs = 0
if nNoiseFreqs == 0
    fp_rate(:) = NaN;
    tn_rate(:) = NaN;
end

% Create results table
resultsTable = struct2table(paramList);
resultsTable.TP_Rate = tp_rate * 100; % In percent
resultsTable.FN_Rate = fn_rate * 100;
resultsTable.FP_Rate = fp_rate * 100;
resultsTable.TN_Rate = tn_rate * 100;

disp('Results Summary Table (Rates in %):');
disp(resultsTable);

% Overall Average Rates (average across parameter sets)
avg_tp_rate = mean(resultsTable.TP_Rate, 'omitnan');
avg_fn_rate = mean(resultsTable.FN_Rate, 'omitnan');
avg_fp_rate = mean(resultsTable.FP_Rate, 'omitnan');
avg_tn_rate = mean(resultsTable.TN_Rate, 'omitnan');

confmat_avg = table(avg_fn_rate, avg_fp_rate, avg_tp_rate, avg_tn_rate, ...
    'VariableNames', {'Avg FN Rate (%)', 'Avg FP Rate (%)', 'Avg TP Rate (%)', 'Avg TN Rate (%)'});
disp('Overall Average Confusion Rates (%):');
disp(confmat_avg);

% Example Plots (Adapt based on what parameters varied)
if numel(param_K_stages) > 1 && numel(param_startWindows) == 1 % Plot vs K
    figure;
    subplot(2, 1, 1);
    plot(resultsTable.K, resultsTable.FP_Rate, 'bo', 'DisplayName', 'FP Rate');
    hold on;
    plot(resultsTable.K, resultsTable.FN_Rate, 'ro', 'DisplayName', 'FN Rate');
    xlabel('Number of Stages (K)');
    ylabel('Rate (%)');
    title(sprintf('False Positive / False Negative Rate vs K (StartWindow=%d)', param_startWindows(1)));
    legend;
    grid on;

    subplot(2, 1, 2);
    plot(resultsTable.K, resultsTable.TP_Rate, 'g^', 'DisplayName', 'TP Rate');
     hold on;
    plot(resultsTable.K, resultsTable.TN_Rate, 'k^', 'DisplayName', 'TN Rate');
    xlabel('Number of Stages (K)');
    ylabel('Rate (%)');
    title(sprintf('True Positive / True Negative Rate vs K (StartWindow=%d)', param_startWindows(1)));
    legend;
    grid on;

elseif numel(param_startWindows) > 1 && numel(param_K_stages) == 1 % Plot vs Start Window
     figure;
    subplot(2, 1, 1);
    plot(resultsTable.startWindow, resultsTable.FP_Rate, 'bo', 'DisplayName', 'FP Rate');
    hold on;
    plot(resultsTable.startWindow, resultsTable.FN_Rate, 'ro', 'DisplayName', 'FN Rate');
    xlabel('Start Window Index');
    ylabel('Rate (%)');
    title(sprintf('False Positive / False Negative Rate vs Start Window (K=%d)', param_K_stages(1)));
    legend;
    grid on;

    subplot(2, 1, 2);
    plot(resultsTable.startWindow, resultsTable.TP_Rate, 'g^', 'DisplayName', 'TP Rate');
     hold on;
    plot(resultsTable.startWindow, resultsTable.TN_Rate, 'k^', 'DisplayName', 'TN Rate');
    xlabel('Start Window Index');
    ylabel('Rate (%)');
    title(sprintf('True Positive / True Negative Rate vs Start Window (K=%d)', param_K_stages(1)));
    legend;
    grid on;
else % Generic box plots if multiple params varied
    figure;
    subplot(2, 2, 1); boxchart(resultsTable.TP_Rate); title('TP Rate (%)'); grid on;
    subplot(2, 2, 2); boxchart(resultsTable.FP_Rate); title('FP Rate (%)'); grid on;
    subplot(2, 2, 3); boxchart(resultsTable.FN_Rate); title('FN Rate (%)'); grid on;
    subplot(2, 2, 4); boxchart(resultsTable.TN_Rate); title('TN Rate (%)'); grid on;
    sgtitle('Distribution of Rates Across All Parameter Sets');
end


% --- Save Results ---
resultsFilename = sprintf('exp_refactored_results_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
try
    save(resultsFilename, 'resultsTable', 'confmat_avg', 'param_startWindows', 'param_K_stages', 'selectedSubjects', 'selectedStimuli');
    fprintf('[%s] Results saved to %s\n', datetime, resultsFilename);
catch ME_save
    warning('Could not save results: %s', ME_save.message);
end

fprintf('[%s] Experiment script finished.\n', datetime);