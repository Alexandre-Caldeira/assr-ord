%% Experiment Runner Script
% Configures, runs, and analyzes experimental data processing.
clearvars; close all; clc;
rng('default'); % For reproducibility if needed

fprintf('===== Starting Experiment Runner =====\n');

%% --- Configuration ---
config = struct();

% Data Loading
config.expDataPath = 'C:\Users\alexa\Desktop\Sinais_EEG\'; % !!! ADJUST AS NEEDED !!!
config.subjects = 1:11;         % Which subjects to process
config.stimuli = 3;             % Which stimuli indices to process (e.g., [3, 4, 5])

% Preprocessing
config.ppc_useZanoteli = true;  % Apply Zanoteli artifact rejection?
config.ppc_useFilter = false;   % Apply Antunes bandpass filter?

% ORD Calculation Parameters (Define parameter space to explore)
config.ord_startWindows = 1;    % Starting window indices
config.ord_K_stages = 2:5;      % Number of stages (K) values
config.ord_method = 'fixedK';   % 'fixedK' or 'fixedStep'
% config.ord_stepSize = 24;     % Only used if method = 'fixedStep'

% ORDTester Parameters
config.tester_desiredAlpha = 0.05;

% Analysis & Plotting
config.plotRateDist = true;
config.plotPareto_FP_TP = true;
config.plotPareto_FN_TP = false; % Example: Add another pareto if desired

% Saving Results
config.saveResults = true;
config.resultsFilename = sprintf('Experiment_Results_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));

%% --- Initialization ---
fprintf('Initializing components...\n');
try
    dtl = DataLoader('exp', 'dataPath', config.expDataPath);
    % Use initial/default Fs, will be updated per file in runExperimentBatch
    ppc = PreProcessor('fs', dtl.fs);
    ordc = ORDCalculator(dtl.fs);
    tester = ORDTester(dtl.signalFrequencies, dtl.noiseFrequencies, ...
                       'desiredAlpha', config.tester_desiredAlpha);
catch ME
    error('Initialization failed: %s', ME.message);
end
fprintf('Initialization complete.\n');

%% --- Run Experiment Batch ---
aggregatedResultsMap = runExperimentBatch(config, dtl, ppc, ordc, tester);

%% --- Post-Processing and Metrics ---
if isempty(aggregatedResultsMap) || aggregatedResultsMap.Count == 0
    error('No results were generated from the experiment batch.');
end

fprintf('Processing aggregated results...\n');
% Convert map to structure format needed by MetricsCalculator
mapKeys = aggregatedResultsMap.keys;
nParamSets = numel(mapKeys);
paramList = repmat(struct(), nParamSets, 1);
aggregatedData = struct(...
    'sumTP', zeros(nParamSets,1), 'sumFN', zeros(nParamSets,1), ...
    'sumFP', zeros(nParamSets,1), 'sumTN', zeros(nParamSets,1), ...
    'totalSignalOpps', zeros(nParamSets,1), 'totalNoiseOpps', zeros(nParamSets,1) ...
    % Add time fields if they were aggregated in runExperimentBatch
);

for i = 1:nParamSets
    key = mapKeys{i};
    data = aggregatedResultsMap(key);
    paramList(i) = data.ordParams; % Store ORD params
    aggregatedData.sumTP(i) = data.sumTP;
    aggregatedData.sumFN(i) = data.sumFN;
    aggregatedData.sumFP(i) = data.sumFP;
    aggregatedData.sumTN(i) = data.sumTN;
    aggregatedData.totalSignalOpps(i) = data.totalSignalOpps;
    aggregatedData.totalNoiseOpps(i) = data.totalNoiseOpps;
    % Assign time data here if aggregated
end

% Create MetricsCalculator instance
metrics = MetricsCalculator(paramList, aggregatedData, ...
                            dtl.nSignalFrequencies, dtl.nNoiseFrequencies, dtl.nChannels); % Use typical values

%% --- Display and Plot Results ---
metrics.displaySummary(config.tester_desiredAlpha); % Display summary and validate FP rate

if config.plotRateDist
    metrics.plotRateDistributions(figure(1));
end
if config.plotPareto_FP_TP
    metrics.plotParetoFront('FP_Rate', 'TP_Rate', figure(2));
end
if config.plotPareto_FN_TP
     metrics.plotParetoFront('FN_Rate', 'TP_Rate', figure(3));
end

%% --- Save Results ---
if config.saveResults
    fprintf('Saving results to %s...\n', config.resultsFilename);
    try
        resultsTable = metrics.getResultsTable();
        confmat_avg = metrics.getAverageResults();
        % Save config along with results
        save(config.resultsFilename, 'resultsTable', 'confmat_avg', 'config');
        fprintf('Results saved successfully.\n');
    catch ME_save
        warning('Could not save results: %s', ME_save.message);
    end
end

fprintf('===== Experiment Runner Finished =====\n');