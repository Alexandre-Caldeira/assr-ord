% --- FILE START: useExample01_SimulationCompare.m ---
% Example script demonstrating simulation, processing, and comparing
% two sequential test strategies (CGST vs. ETS_NDC) using the OO framework.

clear; clc; close all;

fprintf('=== USE EXAMPLE 01: Simulation Mode Comparison ===\n');

% --- Configuration ---
simParams.duration = 120; % Simulation duration in seconds
simParams.fs = 1000;     % Sampling rate
simParams.numChannels = 1; % Number of channels to simulate
% Define simulation noise parameter variations
simNoiseMeans = [-20, -25]; % dB SNR
simNoiseStds = [.1, .1];      % dB Std Dev

% Define MSC analysis parameters to test
mscParamSettings.Method = 'fixedK'; % Use fixed number of stages
mscParamSettings.StartWindows = [1, 5]; % Test starting analysis from window 1 or 5
mscParamSettings.NumStages = [3, 5];   % Test using 5 or 10 stages (K)

% Define Strategy parameters
commonStrategyParams.alpha = 0.05; % Target overall Type I error (for CGST) or per-test alpha (for ETS)
etsStrategyParams = commonStrategyParams;
etsStrategyParams.NDC = 3; % Number of consecutive detections for ETS_NDC

% --- 1. Setup DataManager (Simulation Mode) ---
dataMgr = DataManager('sim'); % Mode 'sim'
dataMgr.Fs = simParams.fs;    % Set sampling rate for simulation
dataMgr.SelectedChannels = 1:simParams.numChannels; % Define channels
% Generate bulk simulation data varying noise mean and std
dataMgr = dataMgr.generateBulkSimulations(simNoiseMeans, simNoiseStds, ...
                                          'duration', simParams.duration, ...
                                          'numChannels', simParams.numChannels);

% --- 2. Setup SignalProcessor ---
processor = SignalProcessor(dataMgr);
processor.ZanoteliEnable = false; % Disable Zanoteli for simple simulation
processor.FilterEnable = false;   % Disable filtering for simple simulation
% Process the generated simulation data (primarily just FFT in this case)
processor = processor.processBulkData(dataMgr);

% --- 3. Setup MSCAnalyzer ---
mscAnalyzer = MSCAnalyzer(mscParamSettings.Method);
% Analyze MSC using the specified parameter ranges
mscAnalyzer = mscAnalyzer.analyzeBulkMSC(processor, ...
                    'StartWindows', mscParamSettings.StartWindows, ...
                    'NumStages', mscParamSettings.NumStages); % Pass ranges

% --- 4 & 5. Setup and Run Tester for CGST Strategy ---
fprintf('\n--- Testing CGST Strategy ---\n');
cgstStrategy = CGST_BetaStrategy();
tester_cgst = SequentialTester(cgstStrategy, dataMgr, mscAnalyzer);
% CGST params only need alpha (M and K come from analysisInfo)
tester_cgst = tester_cgst.runBulkTest(commonStrategyParams); % Use scalar struct
tester_cgst = tester_cgst.calculatePerformance();
tester_cgst = tester_cgst.calculateAverageStoppingTime('Units','stages'); % Calculate avg stop stage

% --- 4 & 5. Setup and Run Tester for ETS_NDC Strategy ---
fprintf('\n--- Testing ETS_NDC Strategy ---\n');
etsStrategy_ndc = ETSStrategy(false); % UseFutility = false
tester_ets_ndc = SequentialTester(etsStrategy_ndc, dataMgr, mscAnalyzer);
% ETS params need alpha and NDC (MM comes from analysisInfo)
tester_ets_ndc = tester_ets_ndc.runBulkTest(etsStrategyParams); % Use scalar struct with NDC
tester_ets_ndc = tester_ets_ndc.calculatePerformance();
tester_ets_ndc = tester_ets_ndc.calculateAverageStoppingTime('Units','stages'); % Calculate avg stop stage

% --- 6. Visualize and Compare Results ---
visualizer_cgst = Visualizer(dataMgr, mscAnalyzer, tester_cgst);
visualizer_ets_ndc = Visualizer(dataMgr, mscAnalyzer, tester_ets_ndc);

fprintf('\n--- Comparison ---\n');
disp('CGST Overall Performance:');
disp(tester_cgst.PerformanceMetrics.Overall);
disp('ETS_NDC Overall Performance:');
disp(tester_ets_ndc.PerformanceMetrics.Overall);

% Plot FPR for both strategies across parameter sets
fig_fpr_cgst = visualizer_cgst.plotFalsePositiveRate();
title(fig_fpr_cgst.CurrentAxes,'FPR vs Parameters (CGST Strategy)');

fig_fpr_ets = visualizer_ets_ndc.plotFalsePositiveRate();
title(fig_fpr_ets.CurrentAxes,'FPR vs Parameters (ETS NDC Strategy)');

% Plot Pareto fronts (using FPR vs Avg Stop Stage)
try
    fig_pareto_cgst = visualizer_cgst.plotParetoFront('XAxisMetric', 'FPR', 'YAxisMetric', 'AverageStoppingTime', 'YAxisUnits', 'stages');
    title(fig_pareto_cgst.CurrentAxes, 'Pareto Front: Avg Stop Stage vs FPR (CGST)');
catch ME_pareto_cgst
    warning('Could not plot CGST Pareto: %s', ME_pareto_cgst.message);
end

try
    fig_pareto_ets = visualizer_ets_ndc.plotParetoFront('XAxisMetric', 'FPR', 'YAxisMetric', 'AverageStoppingTime', 'YAxisUnits', 'stages');
     title(fig_pareto_ets.CurrentAxes, 'Pareto Front: Avg Stop Stage vs FPR (ETS NDC)');
catch ME_pareto_ets
     warning('Could not plot ETS Pareto: %s', ME_pareto_ets.message);
end

% Plot MSC accumulation for a specific case (e.g., MeanIdx 1, StdIdx 1, ParamSet 1,1, Chan 1)
try
    visualizer_cgst.plotMSCandThresholds(1, 1, 1, 1, 1); % Idx1, Idx2, p1Idx, p2Idx, Chan
    title(sprintf('CGST Example Run (MeanIdx=%d, StdIdx=%d, StartWin=%d, K=%d, Chan=%d)', ...
          1,1, mscParamSettings.StartWindows(1), mscParamSettings.NumStages(1), 1));
catch ME_plot
     warning('Could not plot CGST MSC example: %s', ME_plot.message);
end
try
    visualizer_ets_ndc.plotMSCandThresholds(1, 1, 1, 1, 1); % Idx1, Idx2, p1Idx, p2Idx, Chan
    title(sprintf('ETS NDC Example Run (MeanIdx=%d, StdIdx=%d, StartWin=%d, K=%d, Chan=%d)', ...
          1,1, mscParamSettings.StartWindows(1), mscParamSettings.NumStages(1), 1));
catch ME_plot
     warning('Could not plot ETS NDC MSC example: %s', ME_plot.message);
end

fprintf('\n=== USE EXAMPLE 01 COMPLETE ===\n');