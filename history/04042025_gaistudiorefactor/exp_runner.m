%% Experiment Runner Script
% Configures, runs, and analyzes experimental data processing.
clearvars; close all; clc;
rng('default'); % For reproducibility if needed

fprintf('===== Starting Experiment Runner =====\n');

%% --- Configuration ---
config = struct();
config.expDataPath = 'C:\Users\alexa\experimental_data\outros\Sinais_EEG\'; % !!! ADJUST AS NEEDED !!!
config.subjects = 1:11;
config.stimuli = 3;
config.ppc_useZanoteli = true;
config.ppc_useFilter = false;
% --- Choose ONE method ---
config.ord_method = 'fixedK';   % 'fixedK' or 'fixedStep'
config.ord_startWindows = fix(linspace(1,200, 10)); % Only used if method = 'fixedK'
% config.ord_K_stages = 2:5;      % Only used if method = 'fixedK'
% config.ord_method = 'fixedStep';
% config.ord_startWindows = 1; % Typically fixed for fixedStep
% config.ord_stepSize = [24 40 60 100];     % Only used if method = 'fixedStep'
% ---

config.ord_method = 'fixedStep';   % 'fixedK' or 'fixedStep'
config.ord_stepSize = fix(linspace(20,120,50));     % Only used if method = 'fixedStep'

config.tester_desiredAlpha = 0.05;
config.plotBoxcharts = true;
config.plotPareto_FP_TP = true;
config.plotPareto_Time_TP = true; % Plot TP Rate (X) vs AvgDetectTime (Y)
config.paretoColormap = 'copper';
config.plotGroupedPareto = true;  % Plot grouped Pareto (e.g., by noiseMean)
config.groupingVar = 'K_actual'; % Variable to group by
config.addParetoLabels = true;
config.saveResults = true;
config.resultsFilename = sprintf('Experiment_Results_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));

%% --- Initialization ---
fprintf('Initializing components...\n');
try
    dtl = DataLoader('exp', 'dataPath', config.expDataPath);
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
mapKeys = aggregatedResultsMap.keys;
nParamSets = numel(mapKeys);
paramListCell = cell(nParamSets, 1);
aggregatedData = struct(...
    'sumTP', zeros(nParamSets,1), 'sumFN', zeros(nParamSets,1), ...
    'sumFP', zeros(nParamSets,1), 'sumTN', zeros(nParamSets,1), ...
    'totalSignalOpps', zeros(nParamSets,1), 'totalNoiseOpps', zeros(nParamSets,1), ...
    'sumAvgDetTime', zeros(nParamSets,1), 'nAvgDetTime', zeros(nParamSets,1) ... % Add time fields
);
for i = 1:nParamSets
    key = mapKeys{i}; data = aggregatedResultsMap(key);
    ordP = data.ordParams;
    ordP.M = round(mean(data.MList)); % Use average M
    paramListCell{i} = ordP;
    aggregatedData.sumTP(i) = data.sumTP; aggregatedData.sumFN(i) = data.sumFN;
    aggregatedData.sumFP(i) = data.sumFP; aggregatedData.sumTN(i) = data.sumTN;
    aggregatedData.totalSignalOpps(i) = data.totalSignalOpps; aggregatedData.totalNoiseOpps(i) = data.totalNoiseOpps;
    aggregatedData.sumAvgDetTime(i) = data.sumAvgDetTime; % Retrieve aggregated time
    aggregatedData.nAvgDetTime(i) = data.nAvgDetTime;     % Retrieve aggregated count
end
paramList = vertcat(paramListCell{:});

metrics = MetricsCalculator(paramList, aggregatedData, ...
                            dtl.nSignalFrequencies, dtl.nNoiseFrequencies, dtl.nChannels);

%% --- Display and Plot Results ---
metrics.displaySummary(config.tester_desiredAlpha);

metrics.plotAllMetrics( ...
    'includeBoxcharts', config.plotBoxcharts, ...
    'includeParetoFPTP', config.plotPareto_FP_TP, ...
    'includeParetoTimeTP', config.plotPareto_Time_TP, ... % Ensure this is true
    'includeGroupedPareto', config.plotGroupedPareto, ...
    'groupingVar', config.groupingVar, ...
    'paretoColormap', config.paretoColormap, ...
    'addParetoLabels', config.addParetoLabels ...
);

%% --- Save Results ---
if config.saveResults
    fprintf('Saving results to %s...\n', config.resultsFilename);
    try
        resultsTable = metrics.getResultsTable();
        confmat_avg = metrics.getAverageResults();
        save(config.resultsFilename, 'resultsTable', 'confmat_avg', 'config');
        fprintf('Results saved successfully.\n');
    catch ME_save
        warning('Could not save results: %s', ME_save.message);
    end
end

fprintf('===== Experiment Runner Finished =====\n');