% --- FILE START: useExample02_ExperimentalETS.m ---
% Example script demonstrating loading experimental data, processing,
% and applying two ETS strategies (NDC-only vs Pnv) using the OO framework.
% Assumes necessary calibration files (VC_not, parameters) are available.

clear; clc; close all;

fprintf('=== USE EXAMPLE 02: Experimental Data ETS Comparison ===\n');

% --- Configuration ---
% <<< SET THIS PATH to the folder containing Ab50dB.mat etc. >>>
expDataSourcePath = 'C:\Users\alexa\experimental_data\outros\ENTRADAS_PATRICIA\';
% <<< SET THIS PATH to where results should be saved (optional) >>>
resultsSavePath = fullfile(expDataSourcePath, 'OO_Results_Exp');

% --- Check if paths exist ---
if ~exist(expDataSourcePath, 'dir')
    error('Experimental data source path does not exist: %s\nPlease update the path in the script.', expDataSourcePath);
end
if ~exist(resultsSavePath, 'dir')
    fprintf('Creating results directory: %s\n', resultsSavePath);
    mkdir(resultsSavePath);
end

% Select subset of subjects and stimuli for demonstration
subjectsToLoad = [1, 2, 3]; % e.g., 'Ab', 'An', 'Bb'
stimuliToLoad = [3];       % e.g., '50dB' (index 3) -> Corresponds to Mmax=240

% MSC Analysis Parameters (Choose ONE set, matching parameters used for calibration)
mscParamSettings.Method = 'fixedSize';
mscParamSettings.StartWindow = 20; % Example Mmin
mscParamSettings.WindowStepSize = 16;% Example Mstep
Mmax_expected = 240; % Expected Mmax for Stimulus 3 ('50dB')

% Strategy Parameters (MUST correspond to the MSC params and Mmax)
% <<< IMPORTANT: Replace placeholders with actual calibrated values >>>
alpha_per_test = 0.05; % Placeholder
etsParamsNDC.alpha = alpha_per_test;
etsParamsNDC.NDC = 3;    % Placeholder

etsParamsPnv = etsParamsNDC;
% Load VC_not corresponding to Mmax and alpha
% <<< CHECK / UPDATE FILENAME AND PATH >>>
vc_not_filename = sprintf('VC_not_Unitario_Mmax_%d_alfa_%.2f.mat', Mmax_expected, alpha_per_test);
vc_not_filepath = fullfile(expDataSourcePath, vc_not_filename); % Assume it's in the data path
try
    fprintf('Attempting to load VC_not from: %s\n', vc_not_filepath);
    vc_data = load(vc_not_filepath, 'VC_not');
    if ~isvector(vc_data.VC_not) || numel(vc_data.VC_not) < Mmax_expected
         error('Loaded VC_not is not a vector or is too short (length %d, needed >= %d)', numel(vc_data.VC_not), Mmax_expected);
    end
    etsParamsPnv.VC_not = vc_data.VC_not;
    fprintf('Successfully loaded VC_not (length %d).\n', numel(etsParamsPnv.VC_not));
catch ME_vcnot
     warning('Could not load VC_not file "%s". ETS_Pnv strategy will likely fail or use dummy thresholds. Error: %s', vc_not_filename, ME_vcnot.message);
     etsParamsPnv.VC_not = []; % Indicate VC_not is missing
end
% -----------------------------------------------------

% --- 1. Setup DataManager (Experimental Mode) ---
dataMgr = DataManager('exp', expDataSourcePath, resultsSavePath);
dataMgr = dataMgr.loadBulkExperiments(subjectsToLoad, stimuliToLoad);
if dataMgr.Fs ~= 1000; warning('DataManager Fs (%d Hz) differs from expected 1000 Hz.', dataMgr.Fs); end
if dataMgr.ExpSuggestedMMax(stimuliToLoad(1)) ~= Mmax_expected; warning('Configured Mmax_expected (%d) differs from DataManager suggestion (%d) for stimulus %d.', Mmax_expected, dataMgr.ExpSuggestedMMax(stimuliToLoad(1)), stimuliToLoad(1)); end

% --- 2. Setup SignalProcessor ---
processor = SignalProcessor(dataMgr);
processor.ZanoteliEnable = true; % Use Zanoteli preprocessing
processor.FilterEnable = false;  % Disable filtering unless needed
processor = processor.processBulkData(dataMgr);

% --- 3. Setup MSCAnalyzer ---
mscAnalyzer = MSCAnalyzer(mscParamSettings.Method, ...
                          'StartWindow', mscParamSettings.StartWindow, ...
                          'WindowStepSize', mscParamSettings.WindowStepSize);
% Analyze MSC just for this single parameter set
mscAnalyzer = mscAnalyzer.analyzeBulkMSC(processor, ...
                    'StartWindows', mscParamSettings.StartWindow, ...
                    'WindowStepSizes', mscParamSettings.WindowStepSize);

% --- 4 & 5. Setup and Run Tester for ETS_NDC ---
fprintf('\n--- Testing ETS_NDC Strategy ---\n');
etsStrategy_ndc = ETSStrategy(false); % UseFutility = false
tester_ets_ndc = SequentialTester(etsStrategy_ndc, dataMgr, mscAnalyzer);
tester_ets_ndc = tester_ets_ndc.runBulkTest(etsParamsNDC);

% --- 4 & 5. Setup and Run Tester for ETS_Pnv ---
fprintf('\n--- Testing ETS_Pnv Strategy ---\n');
tester_ets_pnv = []; % Initialize
if isempty(etsParamsPnv.VC_not)
    warning('VC_not not loaded. Skipping ETS_Pnv test execution.');
else
    etsStrategy_pnv = ETSStrategy(true); % UseFutility = true
    tester_ets_pnv = SequentialTester(etsStrategy_pnv, dataMgr, mscAnalyzer);
    tester_ets_pnv = tester_ets_pnv.runBulkTest(etsParamsPnv);
end

% --- 6. Visualize Results ---
visualizer_ets_ndc = Visualizer(dataMgr, mscAnalyzer, tester_ets_ndc);
visualizer_ets_pnv = []; % Initialize
if ~isempty(tester_ets_pnv); visualizer_ets_pnv = Visualizer(dataMgr, mscAnalyzer, tester_ets_pnv); end

% Plot MSC vs Thresholds for Subject 1, Stimulus 3 (first loaded), ParamSet (1,1), Channel 1
subjIdxToPlot = subjectsToLoad(1);
stimIdxToPlot = stimuliToLoad(1);
p1IdxToPlot = 1; p2IdxToPlot = 1; % Only one param set analyzed
chanToPlot = 1;
fprintf('\nPlotting example MSC runs for Subject=%d, Stimulus=%d, Channel=%d...\n', subjIdxToPlot, stimIdxToPlot, chanToPlot);

try
    fig_ndc = visualizer_ets_ndc.plotMSCandThresholds(stimIdxToPlot, subjIdxToPlot, p1IdxToPlot, p2IdxToPlot, chanToPlot);
    if ~isempty(fig_ndc); title(fig_ndc.CurrentAxes, sprintf('ETS NDC Run: Subject %d, Stimulus %d, Chan %d', subjIdxToPlot, stimIdxToPlot, chanToPlot)); end
catch ME; warning('Plotting error (ETS NDC): %s', ME.message); end

if ~isempty(visualizer_ets_pnv)
    try
        fig_pnv = visualizer_ets_pnv.plotMSCandThresholds(stimIdxToPlot, subjIdxToPlot, p1IdxToPlot, p2IdxToPlot, chanToPlot);
         if ~isempty(fig_pnv); title(fig_pnv.CurrentAxes, sprintf('ETS Pnv Run: Subject %d, Stimulus %d, Chan %d', subjIdxToPlot, stimIdxToPlot, chanToPlot)); end
    catch ME; warning('Plotting error (ETS Pnv): %s', ME.message); end
else; fprintf('Skipping ETS Pnv example plot as test did not run.\n'); end

% --- Notes on Further Analysis ---
% fprintf('\n--- Experimental Analysis Notes ---\n');
% fprintf('Visualizations show example runs for Subject %d, Stimulus %d, Channel %d.\n', subjIdxToPlot, stimIdxToPlot, chanToPlot);
% fprintf('Performance metrics (TPR, FPR) and confusion matrices require ground truth and are not plotted.\n');
% fprintf('To calculate aggregate TXD, FP, TimeM as in original scripts:\n');
% fprintf('  - Manually define signal/noise frequency bins (e.g., based on "binsM" variable in original .mat files).\n');
% fprintf('  - Access results: decisions = tester_*.TestDecisions{stimIdx, subjIdx, 1, 1}, etc.\n');
% fprintf('  - Aggregate final decisions and stopping stages across subjects.\n');
% fprintf('  - Calculate detection rates and average stopping times.\n');
% 
% fprintf('\n=== USE EXAMPLE 02 COMPLETE ===\n');

% --- 7. Example: Calculate TXD, FP, TimeM (Experimental Data Analysis) ---
fprintf('\n--- Example Experimental Analysis (Aggregated Results) ---\n');

% --- Define Signal/Noise Bins ---
% Attempt to load binsM from one of the original files (assuming consistency)
% If this fails, you MUST define signalBins manually.
signalBins = [];
try
    firstStim = stimuliToLoad(1);
    firstSubj = subjectsToLoad(1);
    subjectID = dataMgr.ExpSubjects{firstSubj};
    stimulusName = dataMgr.ExpStimuli{firstStim};
    fileName = [subjectID, stimulusName, '.mat'];
    filePath = fullfile(dataMgr.DataSourcePath, fileName);
    binsM_data = load(filePath, 'binsM');
    if isfield(binsM_data, 'binsM')
        signalBins = binsM_data.binsM;
        fprintf('Loaded signal bins (binsM) from %s: [%s]\n', fileName, num2str(signalBins));
    else
        warning('Variable "binsM" not found in %s.', fileName);
    end
catch ME_binsM
    warning('Could not load binsM from file %s. Error: %s', fileName, ME_binsM.message);
end

if isempty(signalBins)
    % --- *** MANUALLY DEFINE SIGNAL BINS HERE if loading failed *** ---
    % Example: signalBins = [82, 84, 86, 88, 90, 92, 94, 96]; % Based on SimSignalFrequencies
    error('Signal bins (binsM) could not be loaded and were not defined manually. Stopping analysis.');
end

% Define Noise Bins (Example: +/- N bins around signal, excluding signal)
maxBin = dataMgr.getNumBins();
noiseMargin = 10; % How many bins around signal to consider as potential noise region
noiseBins = [];
for bin_s = signalBins(:)'
    noiseBins = [noiseBins, max(1, bin_s-noiseMargin):min(maxBin, bin_s+noiseMargin)];
end
noiseBins = unique(noiseBins);
noiseBins = setdiff(noiseBins, signalBins); % Exclude signal bins
noiseBins = noiseBins(noiseBins >= 1 & noiseBins <= maxBin); % Ensure valid range
fprintf('Defined noise bins (binsR) around signal bins (+/-%d): %d bins total.\n', noiseMargin, numel(noiseBins));


% --- Aggregate and Display Results for ETS_NDC ---
% fprintf('\nAggregating results for ETS_NDC...\n');
results_ndc = aggregate_exp_results(tester_ets_ndc, signalBins, noiseBins, Mmax_expected);
fprintf('  TXD: %.2f%%\n', results_ndc.TXD);
fprintf('  FP (based on binsR): %.2f%%\n', results_ndc.FP);
fprintf('  Avg. Stop Stage (Signal Bins): %.2f\n', results_ndc.TimeM_Signal);
fprintf('  Avg. Stop Stage (Noise Bins): %.2f\n', results_ndc.TimeM_Noise);

% --- Aggregate and Display Results for ETS_Pnv ---
if ~isempty(tester_ets_pnv)
    fprintf('\nAggregating results for ETS_Pnv...\n');
    results_pnv = aggregate_exp_results(tester_ets_pnv, signalBins, noiseBins, Mmax_expected);
    fprintf('  TXD: %.2f%%\n', results_pnv.TXD);
    fprintf('  FP (based on binsR): %.2f%%\n', results_pnv.FP);
    fprintf('  Avg. Stop Stage (Signal Bins): %.2f\n', results_pnv.TimeM_Signal);
    fprintf('  Avg. Stop Stage (Noise Bins): %.2f\n', results_pnv.TimeM_Noise);
else
    fprintf('\nETS_Pnv test was skipped, no results to aggregate.\n');
end

fprintf('\n--- End Experimental Analysis Example ---\n');

% Helper function (if not already defined or in a separate file)
function out = ifelse(condition, trueVal, falseVal)
    if condition; out = trueVal; else; out = falseVal; end
end

% --- Helper Function to Aggregate Results ---
function results = aggregate_exp_results(tester, signalBins, noiseBins, Mmax_expected)
    results = struct('TXD', NaN, 'FP', NaN, 'TimeM_Signal', NaN, 'TimeM_Noise', NaN);
    if isempty(tester); return; end % Handle skipped tests

    decisionsCell = tester.TestDecisions;
    stoppingCell = tester.StoppingStages;
    if isempty(decisionsCell) || isempty(stoppingCell); return; end

    [nStim, nSubj, nP1, nP2] = size(decisionsCell);
    if nP1 ~= 1 || nP2 ~= 1; warning('Aggregation expects results for only 1 parameter set.'); end

    totalSignalDetections = 0;
    totalSignalOpportunities = 0;
    totalNoiseDetections = 0;
    totalNoiseOpportunities = 0;
    stoppingStagesSignal = [];
    stoppingStagesNoise = [];

    for iStim = 1:nStim % Loop through loaded stimuli (only 1 in this example)
        for jSubj = 1:nSubj % Loop through loaded subjects
            decisions = decisionsCell{iStim, jSubj, 1, 1};
            stopping = stoppingCell{iStim, jSubj, 1, 1};

            if isempty(decisions) || isempty(stopping); continue; end % Skip if test failed

            [nBins, K, nChannels] = size(decisions);
            if K==0 || nChannels==0; continue; end

            final_decision = zeros(nBins, nChannels);
            for chan = 1:nChannels
                for bin = 1:nBins
                    stop_k = stopping(bin, chan);
                    if isnan(stop_k) || stop_k < 1 || stop_k > K; stop_k = K; end
                    stop_k = round(min(stop_k, K));
                    final_decision_val = decisions(bin, stop_k, chan);
                    final_decision(bin, chan) = ifelse(isnan(final_decision_val), 0, final_decision_val);
                end % bin
            end % chan

            % Aggregate counts and stopping stages
            validSigBins = signalBins(signalBins <= nBins);
            validNoiseBins = noiseBins(noiseBins <= nBins);

            if ~isempty(validSigBins)
                decisions_sig = final_decision(validSigBins, :);
                stopping_sig = stopping(validSigBins, :);
                totalSignalDetections = totalSignalDetections + sum(decisions_sig(:) == 1);
                totalSignalOpportunities = totalSignalOpportunities + numel(decisions_sig);
                stoppingStagesSignal = [stoppingStagesSignal; stopping_sig(:)];
            end
            if ~isempty(validNoiseBins)
                decisions_noise = final_decision(validNoiseBins, :);
                stopping_noise = stopping(validNoiseBins, :);
                totalNoiseDetections = totalNoiseDetections + sum(decisions_noise(:) == 1); % FP
                totalNoiseOpportunities = totalNoiseOpportunities + numel(decisions_noise);
                 stoppingStagesNoise = [stoppingStagesNoise; stopping_noise(:)];
            end
        end % subj
    end % stim

    % Calculate Rates and Times
    if totalSignalOpportunities > 0
        results.TXD = 100 * totalSignalDetections / totalSignalOpportunities;
        validStopsSig = stoppingStagesSignal(~isnan(stoppingStagesSignal));
         % Map stages back to window count 'M' using analysisInfo.MM from the strategy
         % This requires getting MM for the *specific* test run, which is complex here.
         % As a fallback, we just average the stage number.
        results.TimeM_Signal = mean(validStopsSig); % Average stopping STAGE for signal
         warning('TimeM_Signal is average STAGE number, not window count (M) or time(s).');
    end
    if totalNoiseOpportunities > 0
        results.FP = 100 * totalNoiseDetections / totalNoiseOpportunities;
         validStopsNoise = stoppingStagesNoise(~isnan(stoppingStagesNoise));
        results.TimeM_Noise = mean(validStopsNoise); % Average stopping STAGE for noise
        warning('TimeM_Noise is average STAGE number, not window count (M) or time(s).');

    end
     % Placeholder for converting stage to time (e.g., if each stage represents Mstep windows of 1s)
     % Need Mstep from mscParamSettings
     % time_per_stage = mscParamSettings.WindowStepSize * 1.0; % Assuming 1s per window
     % results.TimeS_Signal = results.TimeM_Signal * time_per_stage;
end % Helper function
