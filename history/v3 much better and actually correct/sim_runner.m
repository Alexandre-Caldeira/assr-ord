%% Simulation Runner Script
% Configures, runs, and analyzes simulation data processing.
clearvars; close all; clc;
rng('default'); % For reproducible simulations

fprintf('===== Starting Simulation Runner =====\n');

%% --- Configuration ---
config = struct();

% Simulation Parameters (Define parameter space)
config.sim_durations = [180];        % Simulation duration in seconds
config.sim_noiseMeans = [-35 -35];  % Mean SNR dB values
config.sim_noiseStds = [1 5];         % Std Dev SNR dB values

% Preprocessing
config.ppc_useZanoteli = false; % Typically false for sim, but possible
config.ppc_useFilter = false;   % Typically false for sim

% ORD Calculation Parameters (Define parameter space)
config.ord_startWindows = 1:5;  % Starting window indices
config.ord_K_stages = 2:5; %[5 8 10]; %2:5;      % Number of stages (K) values
config.ord_method = 'fixedK';   % 'fixedK' or 'fixedStep'
% config.ord_stepSize = 24;     % Only used if method = 'fixedStep'

% ORDTester Parameters
config.tester_desiredAlpha = 0.05;

% Analysis & Plotting
config.plotRateDist = true;
config.plotTimeDist = true;
config.plotPareto_FP_TP = true;
config.plotPareto_Time_TP = true; % Plot AvgDetectTime vs TP Rate
config.plotPareto_MinTime_TP = false; % Plot AvgMinDetectTime vs TP Rate

% Saving Results
config.saveResults = true;
config.resultsFilename = sprintf('Simulation_Results_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));

%% --- Initialization ---
fprintf('Initializing components...\n');
try
    % Use first sim param values for initial setup
    dtl = DataLoader('sim', 'simDuration', config.sim_durations(1), ...
                     'noiseMean', config.sim_noiseMeans(1), ...
                     'noiseStd', config.sim_noiseStds(1));
    ppc = PreProcessor('fs', dtl.fs); % Fs is fixed for sim
    ordc = ORDCalculator(dtl.fs);
    tester = ORDTester(dtl.signalFrequencies, dtl.noiseFrequencies, ...
                       'desiredAlpha', config.tester_desiredAlpha);
catch ME
    error('Initialization failed: %s', ME.message);
end
fprintf('Initialization complete.\n');

%% --- Run Simulation Batch ---
aggregatedResultsMap = runSimulationBatch(config, dtl, ppc, ordc, tester);

%% --- Post-Processing and Metrics ---
if isempty(aggregatedResultsMap) || aggregatedResultsMap.Count == 0
    error('No results were generated from the simulation batch.');
end

fprintf('Processing aggregated results...\n');

% Convert map to structure format needed by MetricsCalculator
mapKeys = aggregatedResultsMap.keys;
nParamSets = numel(mapKeys);

% Initialize paramList as a CELL array
paramListCell = cell(nParamSets, 1);
paramList = repmat(struct(), nParamSets, 1);
aggregatedData = struct(...
    'sumTP', zeros(nParamSets,1), 'sumFN', zeros(nParamSets,1), ...
    'sumFP', zeros(nParamSets,1), 'sumTN', zeros(nParamSets,1), ...
    'totalSignalOpps', zeros(nParamSets,1), 'totalNoiseOpps', zeros(nParamSets,1), ...
    'sumMinDetTime', zeros(nParamSets,1), 'nMinDetTime', zeros(nParamSets,1), ...
    'sumAvgDetTime', zeros(nParamSets,1), 'nAvgDetTime', zeros(nParamSets,1) ...
);

for i = 1:nParamSets
    key = mapKeys{i};
    data = aggregatedResultsMap(key);

    % Assign to cell array element
    paramListCell{i} = data.ordParams; % Store ORD params struct in cell
    
    aggregatedData.sumTP(i) = data.sumTP;
    aggregatedData.sumFN(i) = data.sumFN;
    aggregatedData.sumFP(i) = data.sumFP;
    aggregatedData.sumTN(i) = data.sumTN;
    aggregatedData.totalSignalOpps(i) = data.totalSignalOpps;
    aggregatedData.totalNoiseOpps(i) = data.totalNoiseOpps;
    aggregatedData.sumMinDetTime(i) = data.sumMinDetTime;
    aggregatedData.nMinDetTime(i) = data.nMinDetTime;
    aggregatedData.sumAvgDetTime(i) = data.sumAvgDetTime;
    aggregatedData.nAvgDetTime(i) = data.nAvgDetTime;
end

% Convert cell array to structure array AFTER the loop
paramList = vertcat(paramListCell{:});

% Create MetricsCalculator instance
metrics = MetricsCalculator(paramList, aggregatedData, ...
                            dtl.nSignalFrequencies, dtl.nNoiseFrequencies, dtl.nChannels);

%% --- Display and Plot Results ---
metrics.displaySummary(config.tester_desiredAlpha); % Display summary and validate FP rate

if config.plotRateDist
    metrics.plotRateDistributions(figure(1));
end
if config.plotTimeDist
    metrics.plotTimeDistributions(figure(2));
end
if config.plotPareto_FP_TP
    metrics.plotParetoFront('FP_Rate', 'TP_Rate', figure(3));
end
if config.plotPareto_Time_TP
    metrics.plotParetoFront('AvgDetectTime_s', 'TP_Rate', figure(4)); % AvgDetectTime is Mean Exam Time
end
if config.plotPareto_MinTime_TP
     metrics.plotParetoFront('AvgMinDetectTime_s', 'TP_Rate', figure(5));
end


%% --- Save Results ---
if config.saveResults
    fprintf('Saving results to %s...\n', config.resultsFilename);
    try
        resultsTable = metrics.getResultsTable();
        confmat_avg = metrics.getAverageResults();
        % Save config along with results
        save(config.resultsFilename, 'resultsTable', 'confmat_avg', 'config');
        fprintf('Results saved successfully.\n', config.resultsFilename);
    catch ME_save
        warning('Could not save results: %s', ME_save.message);
    end
end

fprintf('===== Simulation Runner Finished =====\n');