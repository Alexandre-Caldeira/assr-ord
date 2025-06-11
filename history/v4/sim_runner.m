%% Simulation Runner Script
% Configures, runs, and analyzes simulation data processing.
clearvars; close all; clc;
rng('default'); % For reproducible simulations

fprintf('===== Starting Simulation Runner =====\n');

%% --- Configuration ---
config = struct();

% Simulation Parameters (Define parameter space)
config.sim_durations = [200];        % Simulation duration in seconds
config.sim_noiseMeans = round(linspace(-45,-30,15),0);  % Mean SNR dB values
config.sim_noiseStds = [.1 1] ;%[.1 .5 1 1.5 2 2.5 5].^2;         % Std Dev SNR dB values

% Preprocessing
config.ppc_useZanoteli = false; % Typically false for sim, but possible
config.ppc_useFilter = false;   % Typically false for sim

% ORD Calculation Parameters (Define parameter space)
config.ord_startWindows =  fix(linspace(1,150,5)); % Starting window indices
config.ord_K_stages = 2:5 ;%fix(linspace(2,10,8));      % Number of stages (K) values
config.ord_method = 'fixedK';   % 'fixedK' or 'fixedStep'
% config.ord_stepSize = 24;     % Only used if method = 'fixedStep'

% ORDTester Parameters
config.tester_desiredAlpha = 0.05;

% Analysis & Plotting (Control which plots are generated)
config.plotBoxcharts = true;
config.plotPareto_FP_TP = true;
config.plotPareto_Time_TP = true; % Plot TP Rate (X) vs AvgDetectTime (Y)
config.plotGroupedPareto = true;  % Plot grouped Pareto (e.g., by noiseMean)
config.groupingVar = 'noiseMean'; % Variable to group by
config.paretoColormap = 'abyss';    % Colormap for Pareto plots ('jet', 'parula', 'hsv', etc.)
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
% Returns map where key includes sim params, value has results for that single run
aggregatedResultsMapPerRun = runSimulationBatch(config, dtl, ppc, ordc, tester);

%% --- Post-Processing: Aggregate Across Runs with SAME ORD Params but DIFFERENT Sim Params ---
if isempty(aggregatedResultsMapPerRun) || aggregatedResultsMapPerRun.Count == 0
    error('No results were generated from the simulation batch.');
end

fprintf('Aggregating results across simulation runs...\n');
% This step is crucial for grouping. We create a structure where each row
% represents a unique combination of ORD parameters AND the grouping variable.
finalAggregatedMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
allRunData = aggregatedResultsMapPerRun.values;

for i = 1:numel(allRunData)
    runData = allRunData{i};
    fullParams = runData.fullParams; % Includes M, SW, K, method, duration, noiseMean, noiseStd

    % Create key based on ORD parameters AND the chosen grouping variable
    ordKeyPart = sprintf('SW%d_K%d_Mthd%s', fullParams.startWindow, fullParams.K, fullParams.method);
    groupVal = fullParams.(config.groupingVar); % Get the value of the grouping variable
    finalKey = sprintf('%s_GrpVar%.2f', ordKeyPart, groupVal); % Include group value in key

    % Aggregate into the final map (Summing results for same ORD+Group combination)
    % This shouldn't happen with the current key structure from runSimulationBatch,
    % as each run is unique. This loop essentially reformats the data.
     if finalAggregatedMap.isKey(finalKey)
         % This case should ideally not be hit if runSimulationBatch keys are fully unique.
         % If it is hit, it means multiple runs had identical full parameters.
         aggData = finalAggregatedMap(finalKey);
         aggData.sumTP = aggData.sumTP + runData.sumTP; aggData.sumFN = aggData.sumFN + runData.sumFN;
         aggData.sumFP = aggData.sumFP + runData.sumFP; aggData.sumTN = aggData.sumTN + runData.sumTN;
         aggData.totalSignalOpps = aggData.totalSignalOpps + runData.totalSignalOpps;
         aggData.totalNoiseOpps = aggData.totalNoiseOpps + runData.totalNoiseOpps;
         aggData.sumMinDetTime = aggData.sumMinDetTime + runData.sumMinDetTime;
         aggData.nMinDetTime = aggData.nMinDetTime + runData.nMinDetTime;
         aggData.sumAvgDetTime = aggData.sumAvgDetTime + runData.sumAvgDetTime;
         aggData.nAvgDetTime = aggData.nAvgDetTime + runData.nAvgDetTime;
         finalAggregatedMap(finalKey) = aggData;
     else
         % Create new entry for this ORD+Group combination
         finalAggregatedMap(finalKey) = struct(...
              'params', fullParams,... % Store all parameters for this unique row
              'sumTP', runData.sumTP, 'sumFN', runData.sumFN, 'sumFP', runData.sumFP, 'sumTN', runData.sumTN,...
              'totalSignalOpps', runData.totalSignalOpps, 'totalNoiseOpps', runData.totalNoiseOpps, ...
              'sumMinDetTime', runData.sumMinDetTime, 'nMinDetTime', runData.nMinDetTime, ...
              'sumAvgDetTime', runData.sumAvgDetTime, 'nAvgDetTime', runData.nAvgDetTime ...
              );
     end
end

%% --- Convert Final Aggregated Map to Format for MetricsCalculator ---
fprintf('Converting final aggregated data...\n');
mapKeys = finalAggregatedMap.keys;
nParamSets = numel(mapKeys);
paramListCell = cell(nParamSets, 1);
aggregatedData = struct(...
    'sumTP', zeros(nParamSets,1), 'sumFN', zeros(nParamSets,1), ...
    'sumFP', zeros(nParamSets,1), 'sumTN', zeros(nParamSets,1), ...
    'totalSignalOpps', zeros(nParamSets,1), 'totalNoiseOpps', zeros(nParamSets,1), ...
    'sumMinDetTime', zeros(nParamSets,1), 'nMinDetTime', zeros(nParamSets,1), ...
    'sumAvgDetTime', zeros(nParamSets,1), 'nAvgDetTime', zeros(nParamSets,1) ...
);

for i = 1:nParamSets
    key = mapKeys{i};
    data = finalAggregatedMap(key);
    paramListCell{i} = data.params; % Store the full parameter struct
    aggregatedData.sumTP(i) = data.sumTP; aggregatedData.sumFN(i) = data.sumFN;
    aggregatedData.sumFP(i) = data.sumFP; aggregatedData.sumTN(i) = data.sumTN;
    aggregatedData.totalSignalOpps(i) = data.totalSignalOpps; aggregatedData.totalNoiseOpps(i) = data.totalNoiseOpps;
    aggregatedData.sumMinDetTime(i) = data.sumMinDetTime; aggregatedData.nMinDetTime(i) = data.nMinDetTime;
    aggregatedData.sumAvgDetTime(i) = data.sumAvgDetTime; aggregatedData.nAvgDetTime(i) = data.nAvgDetTime;
end
paramList = vertcat(paramListCell{:}); % Final param list for table

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