% --- FILE START: runCalibrationForExampleParams.m ---

clear; clc; close all;

fprintf('=== Running Calibration for Specific Parameters ===\n');

% --- Configuration (MATCH useExample03) ---
% 1. Simulation Parameters (H0: Noise Only)
nRunsCalib = 10000; % Use a large number for accurate calibration
simParams.duration = 120;
simParams.fs = 1000;
simParams.numChannels = 1;
simNoiseMean = -Inf; simNoiseStd = 0;

% 2. MSC Analysis Parameters (CHOOSE THE *SINGLE* SET TO CALIBRATE)
%    Example: K=10 stages, Target M=12 (needs duration >= 120)
mscParamSettings.Method = 'fixedK';
mscParamSettings.StartWindow = 1;
mscParamSettings.NumStages = 10; % K

% 3. Calibration Targets
nominalAlpha = 0.05; % Target FPR
FP_desejado = nominalAlpha;

% --- Setup DataManager & Processor (Common) ---
dataMgr = DataManager('sim');
dataMgr.Fs = simParams.fs;
dataMgr.SelectedChannels = 1:simParams.numChannels;
dataMgr.SimSignalFrequencies = [];
processor = SignalProcessor(dataMgr);
processor.ZanoteliEnable = false; processor.FilterEnable = false;

% --- Generate H0 MSC Data (as in useExample03) ---
fprintf('Generating %d H0 simulations and calculating MSC...\n', nRunsCalib);
all_msc_data_h0 = cell(nRunsCalib, 1);
actual_K = []; actual_M = []; actual_epochs = [];
for iRun = 1:nRunsCalib
    dataMgr = dataMgr.generateSingleSimulation('duration', simParams.duration, ...
                                               'noiseMean', simNoiseMean, 'noiseStd', simNoiseStd);
    fftData = processor.applyFFT(dataMgr.currentExamTimeSeries);
    mscAnalyzer = MSCAnalyzer(mscParamSettings.Method, ...
                             'StartWindow', mscParamSettings.StartWindow, ...
                             'NumStages', mscParamSettings.NumStages);
    currentParamsMSC = struct('StartWindow', mscParamSettings.StartWindow, ...
                              'NumStages', mscParamSettings.NumStages);
    [mscRun, epochsRun, mRun] = mscAnalyzer.calculateMSC(fftData, currentParamsMSC);
    if isempty(mscRun); continue; end % Skip failed runs
    all_msc_data_h0{iRun} = mscRun;
    if isempty(actual_K); actual_K = size(mscRun, 2); end
    if isempty(actual_M); actual_M = mRun; end
    if isempty(actual_epochs); actual_epochs = epochsRun; end
end
valid_runs = ~cellfun('isempty', all_msc_data_h0);
all_msc_data_h0 = all_msc_data_h0(valid_runs);
nValidRuns = numel(all_msc_data_h0);
if nValidRuns == 0; error('Calibration failed: No valid H0 MSC data generated.'); end
fprintf('Using %d valid H0 runs for calibration.\n', nValidRuns);
fprintf('Achieved K = %d, M = [%s]\n', actual_K, num2str(actual_M));
analysisInfo.K = actual_K; analysisInfo.M = actual_M; analysisInfo.MM = actual_epochs(2:end)-1;

% --- Replicate Calibration Logic (Requires original helper functions) ---
% ASSUMPTION: You have the original .m files for calibration on your path:
%             estimarNDC.m, ETS.m, VC_MSC.m, funcao_custo_v2.m, fmincg.m

% Check for required functions
requiredFiles = {'estimarNDC.m', 'ETS.m', 'VC_MSC.m', 'funcao_custo_v2.m', 'fmincg.m'};
missingFiles = {};
for i=1:length(requiredFiles)
    if isempty(which(requiredFiles{i})); missingFiles{end+1}=requiredFiles{i}; end
end
if ~isempty(missingFiles)
    error('Missing required calibration function files: %s', strjoin(missingFiles, ', '));
end

% Prepare 'ord' matrix expected by original calibration functions
% This requires adapting the MSC data format. Original 'ord' was often (nRuns x Mmax).
% Our data is cell {nRuns} of (nBins x K x nChans).
% We need per-stage MSC values for a representative bin/channel under H0.
% Let's use channel 1 and an arbitrary mid-frequency bin (e.g., bin 100).
targetBin = 100; % Choose a representative noise bin
targetChan = 1;
ord_matrix = zeros(nValidRuns, actual_K);
for iRun = 1:nValidRuns
    mscRun = all_msc_data_h0{iRun};
    if size(mscRun,1) >= targetBin && size(mscRun,3) >= targetChan
       ord_matrix(iRun, :) = squeeze(mscRun(targetBin, :, targetChan));
    else
        warning('Run %d data dimensions [%s] too small for target bin/chan [%d, %d]', iRun, mat2str(size(mscRun)), targetBin, targetChan);
        ord_matrix(iRun, :) = NaN; % Mark as invalid
    end
end
validOrdRows = ~any(isnan(ord_matrix), 2);
ord_matrix = ord_matrix(validOrdRows,:);
fprintf('Prepared "ord" matrix for calibration using %d valid runs, bin %d, chan %d.\n', size(ord_matrix,1), targetBin, targetChan);

% Define the Mmin, Mstep, Mmax equivalent for the chosen parameters
% For fixedK, this mapping is less direct. The original calibration iterates
% through Mmin/Mstep combinations. Here, we have K stages derived from linspace.
% We need analysisInfo.MM which contains the M value for each test point k=1:K.
MM_calib = analysisInfo.MM; % Use the M values at the end of each stage for thresholding
if isempty(MM_calib) || numel(MM_calib)~=actual_K
     error('Cannot determine MM vector needed for calibration functions.');
end

% --- 1. Estimate Optimal NDC (using nominal alpha) ---
fprintf('\nEstimating optimal NDC using nominal alpha=%.3f...\n', nominalAlpha);
try
    % estimarNDC needs Mmin, Mstep, Mmax, but ETS uses MM directly.
    % We might need to adapt estimarNDC or provide dummy Mmin/Mstep/Mmax.
    % Let's assume ETS takes MM directly for now. We need estimarNDC's logic.

    % Simplified logic inspired by estimarNDC: Find NDC yielding FPR <= FP_desejado
    estimatedNDC = NaN;
    maxPossibleNDC = actual_K; % Cannot exceed number of stages
    initialGuessNDC = 1; % Start checking from NDC=1
    minValidNDC = NaN;
    for ndc_test = initialGuessNDC:maxPossibleNDC
        temp_params_ets = struct('alpha', nominalAlpha, 'NDC', ndc_test);
        thresholds_test = etsStrategy_ndc.calculateThresholds(temp_params_ets, analysisInfo); % Use OO class
        % Now apply ETS logic using the *original* ETS function (or replicate its logic)
        dr_ndc = zeros(size(ord_matrix,1), 1);
        for ii = 1:size(ord_matrix,1)
            % ETS(ord_row, MM_vector, alpha_per_test, ndc_value)
            dr_ndc(ii) = ETS(ord_matrix(ii,:)', MM_calib', nominalAlpha, ndc_test); % Pass MM_calib
        end
        currentFPR = mean(dr_ndc);
        fprintf('  NDC=%d -> FPR=%.4f%%\n', ndc_test, currentFPR*100);
        if currentFPR <= FP_desejado
            minValidNDC = ndc_test;
            break; % Found first NDC that meets criteria
        end
    end

    if isnan(minValidNDC)
        warning('Could not find an NDC (up to K=%d) that satisfies FPR <= %.2f%% with nominal alpha. Using NDC=1 for alpha correction.', actual_K, FP_desejado*100);
        estimatedNDC = 1;
    else
        estimatedNDC = minValidNDC; % Use the smallest NDC that worked
        % Optional: Could implement interpolation like original estimarNDC if needed
        fprintf('Estimated optimal NDC = %d\n', estimatedNDC);
    end

catch ME_ndc
     warning('Error during NDC estimation: %s. Using NDC=1 for alpha correction.', ME_ndc.message);
     estimatedNDC = 1;
end

% --- 2. Optimize Alpha (using estimated NDC) ---
fprintf('\nOptimizing alpha using estimated NDC=%d to target FPR=%.2f%%...\n', estimatedNDC, FP_desejado*100);
alfa_teste_optim = nominalAlpha; % Starting point for optimization
options = optimset('MaxIter', 50, 'Display', 'iter'); % Display optimization progress

% Define cost function handle for fmincg
% It needs access to: estimatedNDC, MM_calib, ord_matrix, FP_desejado
costFunc = @(alpha_arg) funcao_custo_v2(alpha_arg, estimatedNDC, MM_calib, ord_matrix, FP_desejado);

try
    % Ensure fmincg is available
    if isempty(which('fmincg'))
        error('fmincg function not found. Ensure it is on the MATLAB path.');
    end
    [calibratedAlpha_ETS_optim, cost] = fmincg(costFunc, alfa_teste_optim, options);

    if isempty(cost) || isnan(calibratedAlpha_ETS_optim)
        warning('Alpha optimization failed. Using nominal alpha.');
        calibratedAlpha_ETS_optim = nominalAlpha;
    else
        fprintf('Optimization complete. Final cost: %g\n', cost(end));
        fprintf('Calibrated Alpha for ETS (NDC=%d) = %.6f\n', estimatedNDC, calibratedAlpha_ETS_optim);
    end
catch ME_optim
     warning('Error during alpha optimization: %s. Using nominal alpha.', ME_optim.message);
     calibratedAlpha_ETS_optim = nominalAlpha;
end

% --- Store Calibrated Values (Example) ---
% In a real scenario, you might save this to a .mat file
% specific to the K, M, NDC parameters.
finalCalibratedAlpha = calibratedAlpha_ETS_optim;
finalEstimatedNDC = estimatedNDC;

fprintf('\n--- Calibration Results ---');
fprintf('\n For K=%d, M=[%s]:', actual_K, num2str(actual_M));
fprintf('\n Estimated NDC = %d', finalEstimatedNDC);
fprintf('\n Calibrated Alpha = %.6f', finalCalibratedAlpha);
fprintf('\n-------------------------\n');

% --- Verification (Optional) ---
% Rerun with calibrated values to check FPR
fprintf('\nVerifying FPR with calibrated alpha=%.6f and NDC=%d...\n', finalCalibratedAlpha, finalEstimatedNDC);
params_verify_ets = struct('alpha', finalCalibratedAlpha, 'NDC', finalEstimatedNDC);
fpr_verified = run_and_get_fpr(etsStrategy_ndc, all_msc_data, params_verify_ets, analysisInfo);
fprintf('Verified FPR: %.3f%%\n', fpr_verified);

fprintf('\n=== CALIBRATION EXAMPLE COMPLETE ===\n');