%% Clear workspace and initialize parameters
clearvars; close all; clc;

% Load object with default exam and reset data
dtl = DataLoader('exp', 'C:\Users\alexa\Desktop\Sinais_EEG\');

dtl.selectedZanoteliSubjects = 1:11; %:numel(dtl.zanoteliSubjects);
dtl.selectedZanoteliStimuli = 3; % 1:numel(dtl.zanoteliStimulusNames)-1; % -1 to remove 'ESP'

% Load data, preprocess and filter
dtl = dtl.loadBulkEEGData();
% ppc = PreProcessor().bulkZanoteliPreprocess(dtl).bulkAntunesFilter(dtl);

% Reset SIGNALS to filtered for display
dtl = dtl.computeBulkFFTs(); % dtl.groupSignals = ppc.groupFilteredSignals;
dtl.age;

%% import parameters from Vaz 2024
caminho = 'C:\Users\alexa\Desktop\Numero_Deteccoes_consecutiva_H\';
vaz_data = load([caminho,'NDC_AlfaCorrigido_Mmax240_alfa_0.05_FPdesejado0.05.mat'], ...
    'alfa_corrigido', 'NDC_minimo','P', 'nRuns');

% vaz_startWindows = vaz_data.P(:,1)';
% vaz_windowSizes = vaz_data.P(:,2)';

% % filter number of stages
vaz_windowSizes = 8:238;
K_stages = zeros(size(vaz_windowSizes));
for idx = 1:numel(K_stages)
    K_stages(idx) = numel(1:vaz_windowSizes(idx):240);
end

% no_go = [find(K_stages >5), find(K_stages<=1)];
% vaz_startWindows(find(K_stages >5)) = [];
% vaz_windowSizes(find(K_stages >5)) = [];
% 
% vaz_windowSizes

vaz_translated_Kstages = flip(unique(K_stages));
vaz_startWindows = 1:240; %240;

%%
ordc = ORDCalculator(dtl).fit_epochs( ...
    stimulusIndices = dtl.selectedZanoteliStimuli,...
    subjectIndices = dtl.selectedZanoteliSubjects,...
    startWindows = vaz_startWindows, ... % windowStepSizes = 52, ...
    K_stages = vaz_translated_Kstages(vaz_translated_Kstages<10),...
    single_or_bulk = 'bulk',...
    lastWindowCalcMethod = 'exactK', ... % maxFromStart, maxFromLast, exactK, fromSizeType
    sizeType = 'fixedSize'... % minToMax, minToFix, withResampling, default = fixedSize
    ... % then, compute on selected channels:
    );


disp('Epochs are computed')
ordc.age()

single_channel = 1;
ordc = ordc.bulk_compute_msc(channels = single_channel);
ordc.age()

%% Test pipeline
% Apply TEST(s) to objective response detector (MSC)
tester = ORDTester(ordc);
tester.desired_alpha = 0.05;

tester = tester.compute_bulk_beta_cgst_decisions();
tester.age()
% disp(tester)

% Comportamento (a implementar)
% O que precisamos testar? Janelas! 
% (qual limite para deteccao, quando parar, como sinalizar)
% 
% Calculamos a posição de inicio e fim das janelas/testes com base em 4 parametros:
% 1. inicio: Instante de inicio do exame
% 2. tamanho: Numero de testes OU Numero de amostras por teste OU Intervalo entre testes;
% 3. paradaT (limita tempo): duração máxima do exame OU numero deamostras/janelas;
% 4. paradaI (limita dados): NDC, futilidade/detecção-CGST, variância das amostras (SNR). 

%% Show results


nParams = numel(tester.epochs(1,1,:));
nChannels = numel(tester.ord_calculator.channels);
tester.FP = zeros(nParams,nChannels);
tester.FN = zeros(nParams,nChannels);
tester.TP = zeros(nParams,nChannels);
tester.TN = zeros(nParams,nChannels);

for stimulusIndex = tester.stimulusIndices
    for subjectIndex = tester.subjectIndices
        for channel_idx = 1:numel(tester.ord_calculator.channels)
            for epoch_param_idx = 1:nParams 
                
                if size(tester.groupFP,1)==1
                    exam_fp = cell2mat(tester.groupFP(1, subjectIndex, epoch_param_idx));
                elseif size(tester.groupFP,2)==1
                    exam_fp = cell2mat(tester.groupFP(stimulusIndex, 1, epoch_param_idx));
                else
                    exam_fp = cell2mat(tester.groupFP(stimulusIndex, subjectIndex, epoch_param_idx));
                end

                if ~isempty(exam_fp)
                    % If there was a detection on any stage, 
                    % assign 1 to that exam, and repeat for all
                    % frequencies. Then, sum all 1s.
                    % + Add previous results from other subjs
                    % To compute rate (pct), divide by #freqs and #subj

                    if size(tester.groupFP,1)==1
                        exam_fn = cell2mat(tester.groupFN(1, subjectIndex, epoch_param_idx));
                        exam_tp = cell2mat(tester.groupTP(1, subjectIndex, epoch_param_idx));
                        exam_tn = cell2mat(tester.groupTN(1, subjectIndex, epoch_param_idx));

                    elseif size(tester.groupFP,2)==1
                        exam_fn = cell2mat(tester.groupFN(stimulusIndex, 1, epoch_param_idx));
                        exam_tp = cell2mat(tester.groupTP(stimulusIndex, 1, epoch_param_idx));
                        exam_tn = cell2mat(tester.groupTN(stimulusIndex, 1, epoch_param_idx));

                    else
                        exam_fn = cell2mat(tester.groupFN(stimulusIndex, subjectIndex, epoch_param_idx));
                        exam_tp = cell2mat(tester.groupTP(stimulusIndex, subjectIndex, epoch_param_idx));
                        exam_tn = cell2mat(tester.groupTN(stimulusIndex, subjectIndex, epoch_param_idx));

                    end
                    
                    

                    tester.FN(epoch_param_idx, channel_idx) = ...
                            sum(any(exam_fn(tester.signalFrequencies,:,channel_idx)>0,2))...
                            + tester.FN(epoch_param_idx, channel_idx);

                    tester.TP(epoch_param_idx, channel_idx) = ...
                            sum(any(exam_tp(tester.signalFrequencies,:,channel_idx)>0,2))...
                            + tester.TP(epoch_param_idx, channel_idx);

                    tester.FP(epoch_param_idx, channel_idx) = ...
                            sum(any(exam_fp(tester.noiseFrequencies,:,channel_idx)>0,2))...
                            + tester.FP(epoch_param_idx, channel_idx);

                    tester.TN(epoch_param_idx, channel_idx) = ...
                            sum(any(exam_tn(tester.noiseFrequencies,:,channel_idx)>0,2))...
                            + tester.TN(epoch_param_idx, channel_idx);
    
                end
            end
        end
    end
end


denom = numel(tester.noiseFrequencies)*numel(tester.subjectIndices)*numel(tester.stimulusIndices);
fp_rate = tester.FP/(denom);

tp_rate = tester.TP/(denom);

fn_rate = tester.FN/(denom);

tn_rate = tester.TN/(denom);

confmat = table( ...
    100*mean(fn_rate(end,:)), ...
    100*mean(fp_rate(end,:)), ...
    100*mean(tp_rate(end,:)), ...
    100*mean(tn_rate(end,:)), ...
    'VariableNames',{'fn','fp','tp','tn'})
