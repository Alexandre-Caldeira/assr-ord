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
% % Author: [Alexandre Gomes Caldeira] 
% % Date: 
% - v0.0: [14 02 2025 10:40] 
%   Project created, exposing incompatibilities...
%
% - v0.1: [07 03 2025 10:45] 
%   DataLoader with FFT functional for exp and sim! :)
%
% - v0.2: [] 
%   Preprocessing functional for exp and sim, signal recompute added
%
% - v0.3: [] 
%   MSC calculator now accepts variable parameters
%
% - v0.4: [] 
%   ORD calculator functional for exp and sim
%
% - v0.5: [] 
%   SHT (single hyp. test [Decisions, Time]) functional for exp and sim 
%
% - v0.6: [] 
%   add Metrics: [Pareto, ConfusionMatrix, Specificity, Repeatability]
%
% - v0.7: [] 
%   GST (group sequential test) functional for exp and sim
%
% - v0.8: [] 
%   GST metrics verified 
%
% - v0.9: []
%
% - v1: [] 
%

%% Clear workspace and initialize parameters
clearvars; 
% close all; 
clc;


% Filter parameters for preprocessing
params.filter.fcLower = 70;                  % Lower cutoff frequency (Hz)
params.filter.fcUpper = params.Fs/2 - 1;     % Upper cutoff frequency (Hz)
params.filter.order   = 8;                   % Butterworth filter order

% Sequential test settings
params.alpha = 0.05;           % Overall false-positive rate
% params.Mmin = []; params.Mmax = [];
params.K     = 7;              % Number of sequential stages
params.M     = params.duration / params.K;  % Window length (seconds) per stage

path_to_params = 'C:\PPGEE\SBEB_CBA_24\CGST_figuras\Sinais_EEG\';
addpath(path_to_params)
load(['NDC_AlfaCorrigido_Mmax' num2str(Mmax) '_alfa_'  num2str(params.alpha) '_FPdesejado' num2str(params.alpha) '.mat'], ...
    'P')
params.P     = P; % parametros = [Min Mstep Mmax alfa_corrigido]

params.testFrequencies = [params.signalFrequencies, params.noiseFrequencies];
params.flagNoise = numel(params.signalFrequencies);

%% Main Pipeline
% dtl = DataLoader('exp');


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
fignum=1;
plotResults(aThresholds, gThresholds, stageMetrics, perf, params,fignum);

% 8. Print summary metrics
fprintf('(FPR,pct.): %.4f\n', perf.FPR);
fprintf('(TPR,pct.): %.4f\n', perf.TPR);
fprintf('(TNR,pct.): %.4f\n', perf.TNR);
fprintf('(FNR,pct.): %.4f\n', perf.FNR);

fprintf('(No decision...,pct.): %.4f\n', sum(decisions==0,'all')/numel(decisions));
fprintf('(Stop: detected.,pct.): %.4f\n', sum(decisions==1,'all')/numel(decisions));
fprintf('(Stop: futile.,pct.): %.4f\n', sum(decisions==-1,'all')/numel(decisions));
