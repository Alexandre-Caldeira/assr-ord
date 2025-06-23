% --- FILE START: useExample03_ThresholdVerification.m ---
% Example script to verify threshold calculation and FPR control
% by comparing results using a nominal alpha vs. a pre-calibrated alpha.

clear; clc; close all;

fprintf('=== USE EXAMPLE 03: Threshold & FPR Verification ===\n');

% --- Configuration ---
% 1. Simulation Parameters (H0: Noise Only)
nRuns = 1; % Number of H0 simulations (Increase for better FPR estimate)
simParams.duration = 120; % Duration resulting in M*K windows (10 stages * 12 wins/stage = 120s)
simParams.fs = 1000;
simParams.numChannels = 1;
simNoiseMean = -Inf; % No signal added (effectively -Inf dB SNR)
simNoiseStd = 0;     % Noise variance doesn't strictly matter for H0 MSC distribution

% 2. MSC Analysis Parameters (CHOOSE ONE SPECIFIC SET)
%    Example: fixedK, K=10 stages, M=12 windows/stage (needs duration >= 120)
mscParamSettings.Method = 'fixedK';
mscParamSettings.StartWindow = 1;
mscParamSettings.NumStages = 10; % K
% To achieve M=12 for K=10 and duration=120: linspace(1, 121, 11) -> diff is 12
% We will verify M later.

% 3. Strategy Parameters
nominalAlpha = 0.05; % The desired overall FPR
etsNDC = 3;          % NDC value for ETS strategy

% 4. <<< IMPORTANT: PRE-CALIBRATED ALPHA VALUES >>>
%    You MUST replace these placeholders with the actual calibrated alpha
%    values obtained from your calibration runs (e.g., funcao_NDC_alfa...)
%    for the *specific* K, M (and potentially NDC) chosen above.
%    If the calibration process yielded different alphas for CGST vs ETS,
%    use the respective values. If it yielded one value based on a specific
%    strategy (e.g., ETS), you might use that for both here for comparison,
%    but be aware CGST's correction might differ.

% Example Placeholder Values (REPLACE THESE!)
calibratedAlpha_CGST = 0.048; % Example: Calibrated alpha for CGST (K=10, M=12)
calibratedAlpha_ETS = 0.045; % Example: Calibrated alpha for ETS (K=10, M=12, NDC=3)

fprintf('--- Configuration ---\n');
fprintf(' H0 Runs: %d\n', nRuns);
fprintf(' Duration: %d s, Fs: %d Hz\n', simParams.duration, simParams.fs);
fprintf(' MSC Params: Method=%s, K=%d (Target), StartWindow=%d\n', ...
        mscParamSettings.Method, mscParamSettings.NumStages, mscParamSettings.StartWindow);
fprintf(' Strategy Params: Nominal Alpha=%.3f, ETS NDC=%d\n', nominalAlpha, etsNDC);
fprintf(' Placeholder Calibrated Alphas: CGST=%.4f, ETS=%.4f (Replace these!)\n', calibratedAlpha_CGST, calibratedAlpha_ETS);
fprintf('--------------------\n');

% --- Setup DataManager ---
dataMgr = DataManager('sim');
dataMgr.Fs = simParams.fs;
dataMgr.SelectedChannels = 1:simParams.numChannels;
dataMgr.SimSignalFrequencies = []; % Ensure NO signal for H0

% --- Generate H0 MSC Data ---
fprintf('Generating %d H0 simulations and calculating MSC...\n', nRuns);
all_msc_data = cell(nRuns, 1);
actual_K = []; % To store the actual K achieved
actual_M = []; % To store the actual M achieved

for iRun = 1:nRuns
    % Generate single noise-only simulation
    dataMgr = dataMgr.generateSingleSimulation('duration', simParams.duration, ...
                                               'noiseMean', simNoiseMean, ...
                                               'noiseStd', simNoiseStd, ...
                                               'signalFreqs', [], ... % Explicitly empty
                                               'numChannels', simParams.numChannels);
    timeSeries = dataMgr.currentExamTimeSeries;

    % Process (just FFT)
    processor = SignalProcessor(dataMgr); % Recreate for correct Fs if needed
    processor.ZanoteliEnable = false;
    processor.FilterEnable = false;
    fftData = processor.applyFFT(timeSeries);

    % Analyze MSC (for the single chosen parameter set)
    mscAnalyzer = MSCAnalyzer(mscParamSettings.Method, ...
                             'StartWindow', mscParamSettings.StartWindow, ...
                             'NumStages', mscParamSettings.NumStages);
    currentParams = struct('StartWindow', mscParamSettings.StartWindow, ...
                           'NumStages', mscParamSettings.NumStages); % Params for calculateMSC
    [mscRun, epochsRun, mRun] = mscAnalyzer.calculateMSC(fftData, currentParams);

    if isempty(mscRun)
        warning('MSC calculation failed for run %d. Skipping.', iRun);
        continue;
    end

    all_msc_data{iRun} = mscRun;

    % Store K and M from the first successful run for validation
    if isempty(actual_K); actual_K = size(mscRun, 2); end
    if isempty(actual_M); actual_M = mRun; end

    if mod(iRun, round(nRuns/10))==0 || iRun==nRuns
         fprintf('  Completed %d / %d runs.\n', iRun, nRuns);
    end
end
fprintf('MSC calculation complete.\n');

% Filter out empty runs
valid_runs = ~cellfun('isempty', all_msc_data);
all_msc_data = all_msc_data(valid_runs);
nValidRuns = numel(all_msc_data);
if nValidRuns < nRuns
    fprintf('Warning: %d runs failed MSC calculation.\n', nRuns - nValidRuns);
    if nValidRuns == 0; error('No valid MSC data generated.'); end
end
fprintf('Using %d valid H0 runs for analysis.\n', nValidRuns);

% --- Verify K and M ---
fprintf('Achieved K = %d stages.\n', actual_K);
if isscalar(actual_M)
    fprintf('Achieved constant M = %d windows per stage.\n', actual_M);
else
    fprintf('Achieved varying M = [%s] windows per stage.\n', num2str(actual_M));
    warning('CGST threshold calculation assumes constant M. Using M(1)=%d.', actual_M(1));
    % Note: If M varies, the calibrated alpha might not be strictly valid for all stages.
end
analysisInfo.K = actual_K;
analysisInfo.M = actual_M; % Pass scalar or vector M
analysisInfo.MM = epochsRun(2:end)-1; % Get MM from the epoch boundaries



% --- Initialize Strategies ---
cgstStrategy = CGST_BetaStrategy();
etsStrategy_ndc = ETSStrategy(false);

% --- Test with NOMINAL Alpha ---
fprintf('\n--- Testing with Nominal Alpha (%.3f) ---\n', nominalAlpha);
params_nominal_cgst = struct('alpha', nominalAlpha);
params_nominal_ets = struct('alpha', nominalAlpha, 'NDC', etsNDC);

fpr_cgst_nominal = run_and_get_fpr(cgstStrategy, all_msc_data, params_nominal_cgst, analysisInfo);
fpr_ets_nominal = run_and_get_fpr(etsStrategy_ndc, all_msc_data, params_nominal_ets, analysisInfo);

% --- Test with CORRECTED Alpha ---
fprintf('\n--- Testing with Calibrated Alpha ---\n');
params_corrected_cgst = struct('alpha', calibratedAlpha_CGST);
params_corrected_ets = struct('alpha', calibratedAlpha_ETS, 'NDC', etsNDC);

fpr_cgst_corrected = run_and_get_fpr(cgstStrategy, all_msc_data, params_corrected_cgst, analysisInfo);
fpr_ets_corrected = run_and_get_fpr(etsStrategy_ndc, all_msc_data, params_corrected_ets, analysisInfo);

% --- Display Comparison ---
fprintf('\n--- FPR Comparison ---\n');
fprintf('Target FPR: %.2f%%\n\n', nominalAlpha*100);
fprintf('Strategy | Nominal Alpha (%.3f) FPR | Calibrated Alpha FPR | Calibrated Alpha Value\n', nominalAlpha);
fprintf('---------|--------------------------|----------------------|------------------------\n');
fprintf('CGST     | %8.3f%%             | %8.3f%%         | %.4f\n', fpr_cgst_nominal, fpr_cgst_corrected, calibratedAlpha_CGST);
fprintf('ETS_NDC  | %8.3f%%             | %8.3f%%         | %.4f\n', fpr_ets_nominal, fpr_ets_corrected, calibratedAlpha_ETS);
fprintf('-----------------------------------------------------------------------------------\n');
if abs(fpr_cgst_corrected - nominalAlpha*100) > 1.0 || abs(fpr_ets_corrected - nominalAlpha*100) > 1.0
    fprintf('NOTE: If corrected FPR is far from target, ensure the provided calibrated alpha values\n');
    fprintf('      *exactly* match the K, M, and NDC parameters used in this script.\n');
end

fprintf('\n=== USE EXAMPLE 03 COMPLETE ===\n');


% --- Helper Function to Run Tests and Calculate FPR ---
function fpr = run_and_get_fpr(strategy, msc_data_cell, strategyParams, analysisInfo)
    n_runs = numel(msc_data_cell);
    fp_count = 0;
    total_tests = 0; % Total bins * channels tested
    fprintf('Calculating thresholds for %s with alpha=%.4f...\n', strategy.Name, strategyParams.alpha);
    try
        thresholds = strategy.calculateThresholds(strategyParams, analysisInfo);
    catch ME
        error('Failed to calculate thresholds: %s', ME.message);
    end
    fprintf('Running %s test on %d runs...\n', strategy.Name, n_runs);
    for r = 1:n_runs
        mscRun = msc_data_cell{r};
        [decisions, ~] = strategy.runTest(mscRun, thresholds);
        if isempty(decisions); continue; end % Skip if test failed

        % Get final decision (1 if detect, 0/-1/NaN otherwise)
        [nBins, K, nChannels] = size(decisions);
        final_decision = decisions(:, K, :); % Decision at the last stage
        % More robust: Use stopping stage if available (assumes stopping stage calc is correct)
        % final_decision = ... see calculatePerformance for logic ...

        fp_count = fp_count + sum(final_decision(:) == 1); % Count detections under H0
        total_tests = total_tests + numel(final_decision);
    end
    if total_tests > 0
        fpr = 100 * fp_count / total_tests;
    else
        fpr = NaN;
    end
    fprintf('  %s: FP Count = %d, Total Tests = %d, FPR = %.3f%%\n', strategy.Name, fp_count, total_tests, fpr);
end