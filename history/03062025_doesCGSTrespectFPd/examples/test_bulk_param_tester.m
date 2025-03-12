%% Clear workspace and initialize parameters
clearvars; close all; clc;

% Load object with default exam and reset data
dtl = DataLoader('exp');

% Define size of vector with randomly selected exams
% and randomly select subjects and stimuli without replacement
% nSubj = 5;
% nStim = 3;
% dtl.selectedZanoteliSubjects = randperm(numel(dtl.zanoteliSubjects),nSubj);
% dtl.selectedZanoteliStimuli = randperm(numel(dtl.zanoteliStimulusNames)-1,nStim); % -1 to remove 'ESP'

% Or choose all
nSubj = 11;
nStim = 5;
dtl.selectedZanoteliSubjects = 1:numel(dtl.zanoteliSubjects);
dtl.selectedZanoteliStimuli = [3 5]; % 1:numel(dtl.zanoteliStimulusNames)-1; % -1 to remove 'ESP'

% Load data, preprocess and filter
dtl = dtl.loadBulkEEGData();
% ppc = PreProcessor().bulkZanoteliPreprocess(dtl);
ppc = PreProcessor().bulkZanoteliPreprocess(dtl).bulkAntunesFilter(dtl);

% Reset SIGNALS to filtered for display
% dtl.groupSignals = ppc.groupProcessedSignals;
dtl.groupSignals = ppc.groupFilteredSignals;
dtl = dtl.computeBulkFFTs();

% M = 40 
ordc = ORDCalculator(dtl).fit_epochs( ...
    stimulusIndices = dtl.selectedZanoteliStimuli,...
    subjectIndices = dtl.selectedZanoteliSubjects,...
    startWindows = [1], ...
    windowStepSizes = 32, ...
    single_or_bulk = 'bulk',...
    lastWindowCalcMethod = 'maxFromStart', ... % maxFromStart, maxFromLast, exactK, fromSizeType
    sizeType = 'fixedSize'... % minToMax, minToFix, withResampling, default = fixedSize
    ... % then, compute on selected channels:
    );

ordc = ordc.bulk_compute_msc(channels = 1);

%% Test pipeline
% Apply TEST(s) to objective response detector (MSC)
% Comportamento:
% O que precisamos testar? Janelas! 
% (qual limite para deteccao, quando parar, como sinalizar)
%
% Calculamos a posição de inicio e fim das janelas/testes com base em 4 parametros:
% 1. inicio: Instante de inicio do exame
% 2. tamanho: Numero de testes OU Numero de amostras por teste OU Intervalo entre testes;
% 3. paradaT (limita tempo): duração máxima do exame OU numero deamostras/janelas;
% 4. paradaI (limita dados): NDC, futilidade/detecção-CGST, variância das amostras (SNR). 
%

tester = ORDTester(ordc).compute_beta_cgst_thresholds();
disp(tester)

%% Show results


