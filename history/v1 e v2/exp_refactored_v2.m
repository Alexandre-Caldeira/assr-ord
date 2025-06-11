%% Clear workspace and initialize parameters
clearvars; close all; clc;
rng('default'); % For reproducible noise frequency generation if needed

fprintf('[%s] Starting Refactored Experiment Script...\n', datetime);

%% --- Configuration ---
% Data Loading
expDataPath = 'C:\Users\alexa\Desktop\Sinais_EEG\'; % !!! ADJUST AS NEEDED !!!
selectedSubjects = 1:11;
selectedStimuli = 3; % Example: 50dB stimulus (index 3)
useAntunesFilter = false; % Set to true to apply filtering
useZanoteliPPC = true;    % Set to true to apply Zanoteli artifact rejection

% ORD Calculation Parameters (Loop Variables)
param_startWindows = 1; % Fixed start window
param_K_stages = 2:5;   % Number of stages (tests) to try
param_method = 'fixedK'; % Use fixed number of stages ('fixedK' or 'fixedStep')
% param_stepSize = 24; % Only used if method = 'fixedStep'

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

%% --- Initialization ---
try
    % Initialize DataLoader with default path, update if needed via Name-Value
    dtl = DataLoader('exp', 'dataPath', expDataPath);
    % Ensure PreProcessor and ORDCalculator use the potentially updated Fs from DataLoader
    % We will update/re-initialize them inside the loop if Fs changes per file,
    % but initialize here to get constants.
    ppc = PreProcessor('fs', dtl.fs);
    ordc = ORDCalculator(dtl.fs);
    tester = ORDTester(dtl.signalFrequencies, dtl.noiseFrequencies, 'desiredAlpha', tester_desiredAlpha);
catch ME
    error('Initialization failed: %s', ME.message);
end

nSignalFreqs = tester.nSignalFrequencies;
nNoiseFreqs = tester.nNoiseFrequencies;
nChannels = dtl.nChannels; % Assuming using all default channels initially

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
        currentFs = dtl.fs; % Get Fs for THIS specific file
        currentChannels = dtl.channels; % Get channels used for THIS file
        currentNChannels = dtl.nChannels; % Get N channels for THIS file

        % Re-initialize/Update dependent objects if Fs changed
        if ppc.fs ~= currentFs
             fprintf('  Updating PreProcessor Fs to %d Hz for this file.\n', currentFs);
             ppc = PreProcessor('fs', currentFs); % Re-init with correct Fs
        end
        if ordc.fs ~= currentFs
             fprintf('  Updating ORDCalculator Fs to %d Hz for this file.\n', currentFs);
             ordc = ORDCalculator(currentFs); % Re-init with correct Fs
        end
        % Tester doesn't depend on Fs directly

        % 2. Preprocess Data
        processedSignals = rawSignals; % Start with raw
        if useZanoteliPPC
            % Use Zanoteli suggested max M, or just use all available after rejection
            maxWindowsForPPC = dtl.zanoteliSuggestedMMax(stim);
             % Ensure maxWindowsForPPC doesn't exceed actual available windows
            maxWindowsForPPC = min(maxWindowsForPPC, size(processedSignals, 2));
            processedSignals = ppc.applyZanoteliPreprocessing(processedSignals, maxWindowsForPPC);
        end
        if useAntunesFilter && ~isempty(processedSignals) % Only filter if data exists
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
                 % Sum outcomes across frequencies AND channels for this configuration
                 tp_matrix = decisionResults.TP(1:nSignalFreqs, :, :); % SigFreqs x K x nChannels
                 fn_matrix = decisionResults.FN(1:nSignalFreqs, :, :);
                 fp_matrix = decisionResults.FP(nSignalFreqs+1:end, :, :); % NoiseFreqs x K x nChannels
                 tn_matrix = decisionResults.TN(nSignalFreqs+1:end, :, :);

                 % Sum occurrences where the outcome happened at ANY stage (k) for each freq/channel pair
                 tp_this_config = sum(any(tp_matrix, 2), [1, 3]); % Sum over dim 1 (freqs) and 3 (channels)
                 fn_this_config = sum(any(fn_matrix, 2), [1, 3]);
                 fp_this_config = sum(any(fp_matrix, 2), [1, 3]);
                 tn_this_config = sum(any(tn_matrix, 2), [1, 3]);

                 % Store aggregated counts for this specific parameter combination
                 paramList = [paramList; currentParams]; %#ok<AGROW>
                 aggregated_TP_counts = [aggregated_TP_counts; tp_this_config]; %#ok<AGROW>
                 aggregated_FN_counts = [aggregated_FN_counts; fn_this_config]; %#ok<AGROW>
                 aggregated_FP_counts = [aggregated_FP_counts; fp_this_config]; %#ok<AGROW>
                 aggregated_TN_counts = [aggregated_TN_counts; tn_this_config]; %#ok<AGROW>

                 % Denominators for rates for THIS configuration
                 nSignalOpportunities = nSignalFreqs * currentNChannels; % Use N channels from THIS file
                 nNoiseOpportunities = nNoiseFreqs * currentNChannels;

                 % Store denominators for later rate calculation per param set
                 aggregated_nSignalOpps = [aggregated_nSignalOpps; nSignalOpportunities]; %#ok<AGROW>
                 aggregated_nNoiseOpps = [aggregated_nNoiseOpps; nNoiseOpportunities]; %#ok<AGROW>

                 % Clear large intermediate data for this param set
                 clear msc decisionResults epochs tp_matrix fn_matrix fp_matrix tn_matrix;

             end % K loop
        end % startWindow loop

        % Clear data for the current subject/stimulus before loading next
        dtl.clearData();
        clear rawSignals processedSignals;
        fprintf('[%s] Finished Subject %d, Stimulus %d.\n', datetime, subj, stim);

    end % Stimulus loop
end % Subject loop

fprintf('\n[%s] All configurations processed. Calculating final metrics...\n', datetime);

%% --- Post-Processing and Visualization ---

if isempty(paramList)
    error('No results were generated. Check parameters and data.');
end

% Prepare aggregated data structure for MetricsCalculator
aggregatedData = struct(...
    'TP', aggregated_TP_counts, ...
    'FN', aggregated_FN_counts, ...
    'FP', aggregated_FP_counts, ...
    'TN', aggregated_TN_counts, ...
    'nSignalOpps', aggregated_nSignalOpps, ... % Pass denominators
    'nNoiseOpps', aggregated_nNoiseOpps ...
);

% Create MetricsCalculator instance
% Pass nSignalFreqs/nNoiseFreqs/nChannels used for the *majority* of tests
% (assuming they don't change drastically, otherwise MetricsCalculator needs adjustment)
metrics = MetricsCalculator(paramList, aggregatedData, nSignalFreqs, nNoiseFreqs, dtl.nChannels); % Use final N channels

% Display Summary
metrics.displaySummary();

% Generate Plots
metrics.plotRateDistributions(figure(1)); % Use figure handles
metrics.plotParetoFront('FP_Rate', 'TP_Rate', figure(2)); % Example Pareto: FP vs TP
% Add more Pareto plots if desired, e.g., FN vs TP
% metrics.plotParetoFront('FN_Rate', 'TP_Rate', figure(3));


% --- Save Results ---
resultsTable = metrics.getResultsTable(); % Get table for saving
confmat_avg = metrics.getAverageResults(); % Get averages for saving
resultsFilename = sprintf('exp_refactored_results_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
try
    % Save relevant config parameters along with results
    save(resultsFilename, 'resultsTable', 'confmat_avg', ...
         'param_startWindows', 'param_K_stages', 'selectedSubjects', 'selectedStimuli', ...
         'nSignalFreqs', 'nNoiseFreqs', 'tester_desiredAlpha', 'useAntunesFilter', 'useZanoteliPPC');
    fprintf('[%s] Results saved to %s\n', datetime, resultsFilename);
catch ME_save
    warning('Could not save results: %s', ME_save.message);
end

fprintf('[%s] Experiment script finished.\n', datetime);