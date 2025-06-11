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
simNoiseMeans = [-10, -15]; % dB SNR
simNoiseStds = [1, .5];      % dB Std Dev

% Define MSC analysis parameters to test
mscParamSettings.Method = 'fixedK'; % Use fixed number of stages
mscParamSettings.StartWindows = [1, 5,10]; % Test starting analysis from window 1 or 5
mscParamSettings.NumStages = [3,4,5,8,10];   % Test using 5 or 10 stages (K)

% Define Strategy parameters
commonStrategyParams.alpha = 0.05; % Target overall Type I error (for CGST) or per-test alpha (for ETS)
etsStrategyParams = commonStrategyParams;
etsStrategyParams.NDC = 10; % Number of consecutive detections for ETS_NDC

% --- 1. Setup DataManager (Simulation Mode) ---
dataMgr = DataManager('sim'); % Mode 'sim'
dataMgr.Fs = simParams.fs;    % Set sampling rate for simulation
dataMgr.SelectedChannels = 1:simParams.numChannels; % Define channels
dataMgr = dataMgr.generateBulkSimulations(simNoiseMeans, simNoiseStds, ...
                                          'duration', simParams.duration, ...
                                          'numChannels', simParams.numChannels);

% --- 2. Setup SignalProcessor ---
processor = SignalProcessor(dataMgr);
processor.ZanoteliEnable = false; % Disable Zanoteli for simple simulation
processor.FilterEnable = false;   % Disable filtering for simple simulation
processor = processor.processBulkData(dataMgr);

% --- 3. Setup MSCAnalyzer ---
mscAnalyzer = MSCAnalyzer(mscParamSettings.Method);
mscAnalyzer = mscAnalyzer.analyzeBulkMSC(processor, ...
                    'StartWindows', mscParamSettings.StartWindows, ...
                    'NumStages', mscParamSettings.NumStages); % Pass ranges

% --- 4 & 5. Setup and Run Tester for CGST Strategy ---
fprintf('\n--- Testing CGST Strategy ---\n');
cgstStrategy = CGST_BetaStrategy();
tester_cgst = SequentialTester(cgstStrategy, dataMgr, mscAnalyzer);
tester_cgst = tester_cgst.runBulkTest(commonStrategyParams); % Use scalar struct
tester_cgst = tester_cgst.calculatePerformance();
tester_cgst = tester_cgst.calculateAverageStoppingTime('Units','stages');

% --- 4 & 5. Setup and Run Tester for ETS_NDC Strategy ---
fprintf('\n--- Testing ETS_NDC Strategy ---\n');
etsStrategy_ndc = ETSStrategy(false); % UseFutility = false
tester_ets_ndc = SequentialTester(etsStrategy_ndc, dataMgr, mscAnalyzer);
tester_ets_ndc = tester_ets_ndc.runBulkTest(etsStrategyParams); % Use scalar struct with NDC
tester_ets_ndc = tester_ets_ndc.calculatePerformance();
tester_ets_ndc = tester_ets_ndc.calculateAverageStoppingTime('Units','stages');

% --- 6. Visualize and Compare Results ---
visualizer_cgst = Visualizer(dataMgr, mscAnalyzer, tester_cgst);
visualizer_ets_ndc = Visualizer(dataMgr, mscAnalyzer, tester_ets_ndc);

fprintf('\n--- Comparison ---\n');
disp('CGST Overall Performance:'); disp(tester_cgst.PerformanceMetrics.Overall);
disp('ETS_NDC Overall Performance:'); disp(tester_ets_ndc.PerformanceMetrics.Overall);

% Plot Overall Performance Bar Chart
try; fig_bar_cgst = visualizer_cgst.plotOverallPerformanceBars(); title(fig_bar_cgst.CurrentAxes, 'Overall Performance Rates (CGST)');
catch ME; warning('Plotting error (CGST Bars): %s', ME.message); end
try; fig_bar_ets = visualizer_ets_ndc.plotOverallPerformanceBars(); title(fig_bar_ets.CurrentAxes, 'Overall Performance Rates (ETS NDC)');
catch ME; warning('Plotting error (ETS Bars): %s', ME.message); end

% Plot Confusion Matrix for first parameter set (P1 Idx = 1, P2 Idx = 1)
p1IdxToPlot = 1; p2IdxToPlot = 1;
fprintf('\nPlotting Confusion Matrix for ParamSet (P1 Idx = %d, P2 Idx = %d)...\n', p1IdxToPlot, p2IdxToPlot);
try; fig_cm_cgst = visualizer_cgst.plotConfusionMatrix(p1IdxToPlot, p2IdxToPlot); title(fig_cm_cgst.CurrentAxes, sprintf('ConfMat (CGST) - %s=%s, %s=%s', mscAnalyzer.analyzedParam1Name, mat2str(mscAnalyzer.analyzedParam1Values(p1IdxToPlot)), mscAnalyzer.analyzedParam2Name, mat2str(mscAnalyzer.analyzedParam2Values(p2IdxToPlot))));
catch ME; warning('Plotting error (CGST CM): %s\nCheck Deep Learning Toolbox.', ME.message); end
try; fig_cm_ets = visualizer_ets_ndc.plotConfusionMatrix(p1IdxToPlot, p2IdxToPlot); title(fig_cm_ets.CurrentAxes, sprintf('ConfMat (ETS NDC) - %s=%s, %s=%s', mscAnalyzer.analyzedParam1Name, mat2str(mscAnalyzer.analyzedParam1Values(p1IdxToPlot)), mscAnalyzer.analyzedParam2Name, mat2str(mscAnalyzer.analyzedParam2Values(p2IdxToPlot))));
catch ME; warning('Plotting error (ETS CM): %s\nCheck Deep Learning Toolbox.', ME.message); end

% Plot FPR across parameter sets
try; fig_fpr_cgst = visualizer_cgst.plotFalsePositiveRate(); title(fig_fpr_cgst.CurrentAxes,'FPR vs Parameters (CGST Strategy)');
catch ME; warning('Plotting error (CGST FPR): %s', ME.message); end
try; fig_fpr_ets = visualizer_ets_ndc.plotFalsePositiveRate(); title(fig_fpr_ets.CurrentAxes,'FPR vs Parameters (ETS NDC Strategy)');
catch ME; warning('Plotting error (ETS FPR): %s', ME.message); end

% Plot Pareto fronts (FPR vs Avg Stop Stage)
try; fig_pareto_cgst = visualizer_cgst.plotParetoFront('XAxisMetric', 'FPR', 'YAxisMetric', 'AverageStoppingTime', 'YAxisUnits', 'stages'); title(fig_pareto_cgst.CurrentAxes, 'Pareto Front: Avg Stop Stage vs FPR (CGST)');
catch ME; warning('Plotting error (CGST Pareto): %s', ME.message); end
try; fig_pareto_ets = visualizer_ets_ndc.plotParetoFront('XAxisMetric', 'FPR', 'YAxisMetric', 'AverageStoppingTime', 'YAxisUnits', 'stages'); title(fig_pareto_ets.CurrentAxes, 'Pareto Front: Avg Stop Stage vs FPR (ETS NDC)');
catch ME; warning('Plotting error (ETS Pareto): %s', ME.message); end

% Plot MSC example run
fprintf('\nPlotting example MSC runs for MeanIdx=1, StdIdx=1, P1 Idx = %d, P2 Idx = %d, Chan=1...\n', p1IdxToPlot, p2IdxToPlot);
try; visualizer_cgst.plotMSCandThresholds(1, 1, p1IdxToPlot, p2IdxToPlot, 1); title(sprintf('CGST Example Run (MIdx=%d,SIdx=%d,%s=%s,%s=%s,Ch=1)', 1,1, mscAnalyzer.analyzedParam1Name, mat2str(mscAnalyzer.analyzedParam1Values(p1IdxToPlot)), mscAnalyzer.analyzedParam2Name, mat2str(mscAnalyzer.analyzedParam2Values(p2IdxToPlot))));
catch ME; warning('Plotting error (CGST MSC): %s', ME.message); end
try; visualizer_ets_ndc.plotMSCandThresholds(1, 1, p1IdxToPlot, p2IdxToPlot, 1); title(sprintf('ETS NDC Example Run (MIdx=%d,SIdx=%d,%s=%s,%s=%s,Ch=1)', 1,1, mscAnalyzer.analyzedParam1Name, mat2str(mscAnalyzer.analyzedParam1Values(p1IdxToPlot)), mscAnalyzer.analyzedParam2Name, mat2str(mscAnalyzer.analyzedParam2Values(p2IdxToPlot))));
catch ME; warning('Plotting error (ETS MSC): %s', ME.message); end

fprintf('\n=== USE EXAMPLE 01 COMPLETE ===\n');