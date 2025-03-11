%% Clear workspace and initialize parameters
clearvars; close all; clc;

% Load object with default exam and reset data
dtl = DataLoader('exp');

%% Test pipeline
% Define size of vector with randomly selected exams
% and randomly select subjects and stimuli without replacement
% nSubj = 5;
% nStim = 3;
% dtl.selectedZanoteliSubjects = randperm(numel(dtl.zanoteliSubjects),nSubj);
% dtl.selectedZanoteliStimuli = randperm(numel(dtl.zanoteliStimulusNames),nStim);

% Or choose all
% nSubj = 11;
% nStim = 5;
% dtl.selectedZanoteliSubjects = 1:numel(dtl.zanoteliSubjects);
% dtl.selectedZanoteliStimuli = 1:numel(dtl.zanoteliStimulusNames)-1; % -1 to remove 'ESP'

% Load data, preprocess and filter
dtl = dtl.loadBulkEEGData();
ppc = PreProcessor().bulkZanoteliPreprocess(dtl).bulkAntunesFilter(dtl);

% Reset SIGNALS to filtered for display
dtl.groupSignals = ppc.groupFilteredSignals;
dtl = dtl.computeBulkFFTs();

% Compute objective response detector (MSC)
% Comportamento:
% O que precisamos definir? Janelas! 
% (Onde cada começa, quantas amostras têm.)
%
% Calculamos a posição de inicio e fim delas com base em 4 parametros:
% 1. inicio: Instante de inicio do exame
% 2. tamanho: Numero de testes OU Numero de amostras por teste OU Intervalo entre testes;
% 3. paradaT (limita tempo): duração máxima do exame OU numero deamostras/janelas;
% 4. paradaI (limita dados): NDC, futilidade/detecção-CGST, variância das amostras (SNR). 
%
% 

% compute_bulk_msc( ...
% startWindows = 1:5:21, ...
% windowStepSizes = [5 18 24 32], ...
% channels = [1 8 16], ...
% lastWindowCalcMethod = 'flexible', ... % maxFromStart, maxFromLast, exactK, flexible
% stepType = 'fixedSize'... % minToK, minToMax, minToFix, withResampling, default = fixedSize
%     ...
%     );

% Specify epoch parameters and compute MSCs accordingly
ordc = ORDCalculator(dtl).fit_epochs( ...
    startWindows = [1 5 15], ...
    windowStepSizes = [5 18 24 32], ...
    channels = [1 8 16], ...
    lastWindowCalcMethod = 'maxFromStart', ... % maxFromStart, maxFromLast, exactK, fromSizeType
    sizeType = 'fixedSize'... % minToMax, minToFix, withResampling, default = fixedSize
    ).compute_bulk_msc();

%% Show results
% Show object
dtl.age()
disp(ordc)

filtered_freq_ord = ordc.latestMSC;

lead_name = dtl.zanoteliLeads(random_electrode);
random_epoch = random_epoch+2; % add 2 seconds removed during preprocessing

exam = figure(3);
xlabel('Frequency [Hz]')
ylabel('Voltage [?]')
xlim([0,dtl.nBins])
ylim([0,min(1.10*max(ordc.MSC(:,:,random_electrode),[],'all'),1.09)])
grid on
hold on
set(exam, 'WindowState', 'maximized');
nWindows = size(ordc.MSC,2);
for window_index = 1:nWindows
    title_str =['Window #', num2str(window_index),' of ',...
        num2str(nWindows),': MSC filtered from ' cell2mat(lead_name)];
    title(title_str)

    s = stem(squeeze(ordc.MSC(:,window_index,random_electrode)), ...
        'filled', 'MarkerSize',3,'LineWidth',0.1, 'Color','#0072BD');

    t = stem(dtl.signalFrequencies, ...
        squeeze(ordc.MSC(dtl.signalFrequencies,window_index,random_electrode)), ...
        'filled', 'MarkerSize',3,'LineWidth',0.1, 'Color','red');

    drawnow
    if window_index<nWindows
        pause(1.3)
        delete(s)
        delete(t)
    end
end

    

