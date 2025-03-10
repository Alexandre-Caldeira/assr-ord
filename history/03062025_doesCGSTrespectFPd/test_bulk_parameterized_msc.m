%% Clear workspace and initialize parameters
clearvars; 
close all; 
clc;

%% Test pipeline
% Load object with default exam and reset data
dtl = DataLoader('exp');

% Define size of vector with randomly selected exams
nSubj = 7;
nStim = 4;

% Randomly select subjects and stimuli without replacement
dtl.selectedZanoteliSubjects = randperm(numel(dtl.zanoteliSubjects),nSubj);
dtl.selectedZanoteliStimuli = randperm(numel(dtl.zanoteliStimulusNames),nStim);

% Or choose all
% dtl.selectedZanoteliSubjects = 1:numel(dtl.zanoteliSubjects);
% dtl.selectedZanoteliStimuli = 1:numel(dtl.zanoteliStimulusNames);

% Load data and compute FFT
dtl = dtl.loadBulkEEGData().computeBulkFFTs();

% Preprocess and filter all data
ppc = ppc.bulkZanoteliPreprocess().bulkAntunesFilter();

% Compute objective response detector (MSC)
ordc = ORDCalculator(dtl);
ordc = ordc.compute_msc();
