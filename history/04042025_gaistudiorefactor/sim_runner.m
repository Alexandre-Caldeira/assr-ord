%% Simulation Runner Script
% Configures, runs, and analyzes simulation data processing.
clearvars; close all; clc;
rng('default'); % For reproducible simulations

fprintf('===== Starting Simulation Runner =====\n');

%% --- Configuration ---
config = struct();

% Simulation Parameters (Define parameter space)
config.sim_durations = [480];        % Simulation duration in seconds
config.sim_noiseMeans = round(linspace(-20,-25,-15),0);  % Mean SNR dB values
config.sim_noiseStds = [.5 .5] ;%[.1 .5 1 1.5 2 2.5 5].^2;         % Std Dev SNR dB values

% Preprocessing
config.ppc_useZanoteli = false; % Typically false for sim, but possible
config.ppc_useFilter = false;   % Typically false for sim

% ORD Calculation Parameters (Define parameter space)
% config.ord_startWindows =  fix(linspace(1,150,5)); % Starting window indices
config.ord_K_stages = [2 4 5]; %fix(linspace(2,10,8));      % Number of stages (K) values
config.ord_method = 'fixedStep';   % 'fixedK' or 'fixedStep'
config.ord_stepSize = [24 40 60 100];     % Only used if method = 'fixedStep'

% ORDTester Parameters
config.tester_desiredAlpha = 0.05;

% Analysis & Plotting (Control which plots are generated)
config.plotBoxcharts = true;
config.plotPareto_FP_TP = true;
config.plotPareto_Time_TP = true; % Plot TP Rate (X) vs AvgDetectTime (Y)
config.plotGroupedPareto = true;  % Plot grouped Pareto (e.g., by noiseMean)
config.groupingVar = 'K_actual'; % Variable to group by
config.paretoColormap = 'copper';    % Colormap for Pareto plots ('jet', 'parula', 'hsv', 'turbo', etc.)
config.addParetoLabels = true;    % Add labels to Pareto points?

% Saving Results
config.saveResults = true;
config.resultsFilename = sprintf('Simulation_Results_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));

%% --- Initialization ---
fprintf('Initializing components...\n');
try
    dtl = DataLoader('sim', 'simDuration', config.sim_durations(1), ...
                     'noiseMean', config.sim_noiseMeans(1), ...
                     'noiseStd', config.sim_noiseStds(1));
    ppc = PreProcessor('fs', dtl.fs);
    ordc = ORDCalculator(dtl.fs);
    tester = ORDTester(dtl.signalFrequencies, dtl.noiseFrequencies, ...
                       'desiredAlpha', config.tester_desiredAlpha);
catch ME
    error('Initialization failed: %s', ME.message);
end
fprintf('Initialization complete.\n');

%% --- Run Simulation Batch ---
aggregatedResultsMapPerRun = runSimulationBatch(config, dtl, ppc, ordc, tester);

%% --- Post-Processing: Create Table for MetricsCalculator ---
if isempty(aggregatedResultsMapPerRun) || aggregatedResultsMapPerRun.Count == 0
    error('No results were generated from the simulation batch.');
end

fprintf('Converting results map to table format...\n');
allRunData = aggregatedResultsMapPerRun.values;
nParamSets = numel(allRunData);
paramListCell = cell(nParamSets, 1);
aggregatedData = struct(...
    'sumTP', zeros(nParamSets,1), 'sumFN', zeros(nParamSets,1), ...
    'sumFP', zeros(nParamSets,1), 'sumTN', zeros(nParamSets,1), ...
    'totalSignalOpps', zeros(nParamSets,1), 'totalNoiseOpps', zeros(nParamSets,1), ...
    'sumMinDetTime', zeros(nParamSets,1), 'nMinDetTime', zeros(nParamSets,1), ...
    'sumAvgDetTime', zeros(nParamSets,1), 'nAvgDetTime', zeros(nParamSets,1) ...
);

for i = 1:nParamSets
    runData = allRunData{i};
    paramListCell{i} = runData.fullParams; % Store the full parameter struct
    aggregatedData.sumTP(i) = runData.sumTP; aggregatedData.sumFN(i) = runData.sumFN;
    aggregatedData.sumFP(i) = runData.sumFP; aggregatedData.sumTN(i) = runData.sumTN;
    aggregatedData.totalSignalOpps(i) = runData.totalSignalOpps; aggregatedData.totalNoiseOpps(i) = runData.totalNoiseOpps;
    aggregatedData.sumMinDetTime(i) = runData.sumMinDetTime; aggregatedData.nMinDetTime(i) = runData.nMinDetTime;
    aggregatedData.sumAvgDetTime(i) = runData.sumAvgDetTime; aggregatedData.nAvgDetTime(i) = runData.nAvgDetTime;
end
paramList = vertcat(paramListCell{:});

%% --- Create MetricsCalculator Instance ---
metrics = MetricsCalculator(paramList, aggregatedData, ...
                            dtl.nSignalFrequencies, dtl.nNoiseFrequencies, dtl.nChannels);

%% --- Display and Plot Results ---
metrics.displaySummary(config.tester_desiredAlpha);

metrics.plotAllMetrics( ...
    'includeBoxcharts', config.plotBoxcharts, ...
    'includeParetoFPTP', config.plotPareto_FP_TP, ...
    'includeParetoTimeTP', config.plotPareto_Time_TP, ...
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
        fprintf('Results saved successfully.\n', config.resultsFilename);
    catch ME_save
        warning('Could not save results: %s', ME_save.message);
    end
end

fprintf('===== Simulation Runner Finished =====\n');