%% Experiment Runner Script
% Configures, runs, and analyzes experimental data processing.
clearvars; close all; clc;
rng('default'); % For reproducibility if needed

fprintf('===== Starting Experiment Runner =====\n');

%% --- Configuration ---
config = struct();
config.expDataPath = 'C:\PPGEE\SBEB_CBA_24\CGST_figuras\Sinais_EEG'; % !!! ADJUST AS NEEDED !!!
config.subjects = 1:11;
config.stimuli = 3;
config.ppc_useZanoteli = true;
config.ppc_useFilter = false;
config.ord_startWindows = 1:5;
config.ord_K_stages = 2:5;
config.ord_method = 'fixedK';
config.tester_desiredAlpha = 0.05;
config.plotBoxcharts = true; % Use boxcharts
config.plotPareto_FP_TP = true;
config.plotPareto_Time_TP = true; % Plot TP Rate (X) vs AvgDetectTime (Y)
config.plotPareto_FN_TP = false;
config.paretoColormap = 'abyss';
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
% Returns map where key is ORD params, value aggregates over subj/stim
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
    'totalSignalOpps', zeros(nParamSets,1), 'totalNoiseOpps', zeros(nParamSets,1) ...
);
for i = 1:nParamSets
    key = mapKeys{i}; data = aggregatedResultsMap(key);
    % Calculate average M for this ORD setting and add to params
    ordP = data.ordParams;
    ordP.M = round(mean(data.MList)); % Use average M
    paramListCell{i} = ordP;
    aggregatedData.sumTP(i) = data.sumTP; aggregatedData.sumFN(i) = data.sumFN;
    aggregatedData.sumFP(i) = data.sumFP; aggregatedData.sumTN(i) = data.sumTN;
    aggregatedData.totalSignalOpps(i) = data.totalSignalOpps; aggregatedData.totalNoiseOpps(i) = data.totalNoiseOpps;
end
paramList = vertcat(paramListCell{:});

metrics = MetricsCalculator(paramList, aggregatedData, ...
                            dtl.nSignalFrequencies, dtl.nNoiseFrequencies, dtl.nChannels);

%% --- Display and Plot Results ---
metrics.displaySummary(config.tester_desiredAlpha);

metrics.plotAllMetrics( ...
    'includeBoxcharts', config.plotBoxcharts, ...
    'includeParetoFPTP', config.plotPareto_FP_TP, ...
    'includeParetoTimeTP', true, ... % No time data
    'includeGroupedPareto', false, ... % No grouping
    'paretoColormap', config.paretoColormap, ...
    'addParetoLabels', config.addParetoLabels ...
    'groupingVar', '', ...
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