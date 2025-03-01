% function SinalEEG_Analysis() 
% SinalEEG_Analysis 
% % This unified MATLAB application: 
% • Loads or simulates EEG data 
% • Preprocesses the signal (DC removal, bandpass filtering) 
% • Computes FFT and magnitude-squared coherence (MSC) 
% • Computes beta-distribution based sequential test thresholds 
% • Runs a multi-stage sequential test with early stopping 
% • Evaluates performance metrics and visualizes the results 
% % Developed based on modern software engineering practices. 
% % Author: [Alexandre Gomes Caldeira] % Date: [14 02 2025 10:40]

%% Clear workspace and initialize parameters
clearvars; 
% close all; 
clc;

% Create a central parameter structure
params = struct();
params.mode = 'simulation';  % Options: 'simulation' or 'real'

% Sampling and FFT parameters
params.Fs   = 1000;                  % Sampling frequency (Hz)
params.nfft = params.Fs;             % Number of FFT points (1 sec window)
params.nBins = floor(params.Fs/2)+1; % ...

% Filter parameters for preprocessing
params.filter.fcLower = 70;                  % Lower cutoff frequency (Hz)
params.filter.fcUpper = params.Fs/2 - 1;     % Upper cutoff frequency (Hz)
params.filter.order   = 8;                   % Butterworth filter order

% Simulation settings (used when mode=='simulation')
% params.sim.NumTests    = 1e5;    % Number of Monte Carlo tests
params.sim.signalSNR   = -5;     % Signal-to-noise ratio in dB
params.duration = 120; % seconds
params.nChannels = 3;

% Sequential test settings
params.alpha = 0.10;           % Overall false-positive rate
% params.Mmin = []; params.Mmax = [];
params.K     = 7;              % Number of sequential stages
params.M     = params.duration / params.K;  % Window length (seconds) per stage
params.signalFrequencies = [82    84    86    88    90    92    94    96]; 
params.noiseFrequencies  = 300:407;

params.testFrequencies = [params.signalFrequencies, params.noiseFrequencies];
params.flagNoise = numel(params.signalFrequencies);

% Data file path (if using real EEG data)
% e.g., 'C:\Data\EEG\subject01.mat'
params.data.path = 'C:\PPGEE\SBEB_CBA_24\CGST_figuras\Sinais_EEG\Ab30dB.mat';  

%% Main Pipeline

% 1. Data Loading & Preprocessing
if strcmpi(params.mode, 'simulation')
    fprintf('Running simulation mode...\n');
    [signals, params] = genSimulatedSignals(params);
else
    fprintf('Loading real EEG data...\n');
    signals = loadEEGData(params.data.path);
    signals = preprocessSignal(signals, params);
end

% 2. Compute FFT of the signals
fftSignals = computeFFT(signals, params);

% 3. Compute MSC (here a simplified average per stage)
MSCvalues = computeMSC(fftSignals, params);

% 4. Compute sequential test thresholds
[aThresholds, gThresholds] = computeBetaThresholds(params.K, params.M, params);
params.aThresholds = aThresholds;
params.gThresholds = gThresholds;

% 5. Run the sequential test (decision process)
[decisions, stageMetrics] = sequentialTest(MSCvalues, params);

% 6. Evaluate performance metrics (e.g., FPR and TPR)
% Comment: 
% gets Attempt to grow array along ambiguous dimension.
% when numel(noiseFrequencies) > numel(signalFrequencies)
% but not otherwise
% why?
perf = performanceMetrics(decisions, params);

% 7. Plot the results
% fignum=2;
% plotResults(aThresholds, gThresholds, stageMetrics, perf, params,fignum);

% 8. Print summary metrics
fprintf('(FPR,pct.): %.4f\n', perf.FPR);
fprintf('(TPR,pct.): %.4f\n', perf.TPR);
fprintf('(TNR,pct.): %.4f\n', perf.TNR);
fprintf('(FNR,pct.): %.4f\n', perf.FNR);

fprintf('(No decision...,pct.): %.4f\n', sum(isnan(decisions),'all')/numel(decisions));
fprintf('(Stop: detected.,pct.): %.4f\n', sum(decisions==1,'all')/numel(decisions));
fprintf('(Stop: futile.,pct.): %.4f\n', sum(decisions==-1,'all')/numel(decisions));

% end

%% Data Loading and Preprocessing Functions
% 
% function signals = loadEEGData(filepath) % loadEEGData loads EEG data from the specified file. % The file should contain a variable named 'x' representing the EEG data. if isempty(filepath) error('EEG data file path not specified.'); end data = load(filepath); if isfield(data, 'x') signals = data.x; else error('The EEG data file must contain variable "x".'); end end
% 
% function procSignals = preprocessSignal(signals, params) % preprocessSignal removes the DC offset and applies a Butterworth bandpass filter. procSignals = signals - mean(signals, 2); [b, a] = butter(params.filter.order, ... [params.filter.fcLower, params.filter.fcUpper] / (params.Fs/2)); procSignals = filtfilt(b, a, procSignals')'; end
% 
% function fftSignals = computeFFT(signals, nfft) % computeFFT computes the FFT along each window (row-wise) and returns only the positive frequencies. fftSignals = fft(signals, nfft, 2); fftSignals = fftSignals(:, 1:floor(nfft/2)+1); end
% 
%% Feature Extraction Functions
% 
% function MSCvalues = computeMSC(fftSignals, K) % computeMSC calculates a simplified magnitude-squared coherence (MSC) value per stage. totalBins = size(fftSignals, 2); windowsPerStage = floor(totalBins / K); MSCvalues = zeros(1, K); for k = 1:K idxStart = (k-1)*windowsPerStage + 1; idxEnd = idxStart + windowsPerStage - 1; stageData = abs(fftSignals(:, idxStart:idxEnd)).^2; MSCvalues(k) = mean(stageData(:)); end end
% 
%% Threshold and Sequential Test Functions
% 
% function [aThresholds, gThresholds] = computeBetaThresholds(K, M, alpha) % computeBetaThresholds calculates sequential detection and futility thresholds. aThresholds = zeros(1, K); gThresholds = zeros(1, K);