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
%%
% M = 40 
K_stages = 3;
ordc = ORDCalculator(dtl).fit_epochs( ...
    stimulusIndices = dtl.selectedZanoteliStimuli,...
    subjectIndices = dtl.selectedZanoteliSubjects,...
    startWindows = [1 5], ... % windowStepSizes = 52, ...
    K_stages = K_stages,...
    single_or_bulk = 'bulk',...
    lastWindowCalcMethod = 'exactK', ... % maxFromStart, maxFromLast, exactK, fromSizeType
    sizeType = 'fixedSize'... % minToMax, minToFix, withResampling, default = fixedSize
    ... % then, compute on selected channels:
    );

single_channel = 1;
ordc = ordc.bulk_compute_msc(channels = single_channel);
dtl.age()

%% Test pipeline
% Apply TEST(s) to objective response detector (MSC)
tester = ORDTester(ordc);
tester = tester.compute_bulk_beta_cgst_decisions();


%%
% tester = tester.validateDetectionThresholds(dtl);
disp(tester)

% Comportamento:
% O que precisamos testar? Janelas! 
% (qual limite para deteccao, quando parar, como sinalizar)
%
% Calculamos a posição de inicio e fim das janelas/testes com base em 4 parametros:
% 1. inicio: Instante de inicio do exame
% 2. tamanho: Numero de testes OU Numero de amostras por teste OU Intervalo entre testes;
% 3. paradaT (limita tempo): duração máxima do exame OU numero deamostras/janelas;
% 4. paradaI (limita dados): NDC, futilidade/detecção-CGST, variância das amostras (SNR). 

%% Show results

check = tester.FP + tester.FN + tester.TP + tester.TN;
check_tested = check(tester.allTestFrequencies,:,:);
check_untested = check;
check_untested(tester.allTestFrequencies,:,:) = [];

% decisao_por_exame = zeros(params.dataloader.noiseFrequencies,size(tester.FP,3));
% fp_por_exame = squeeze(any(tester.FP(tester.dataloader.noiseFrequencies,:,:)>0,2));
% fp_medio = 100*mean(fp_por_exame ,'all')
% 
% tn_por_exame = squeeze(any(tester.FP(tester.dataloader.noiseFrequencies,:,:)>0,2));
% tn_medio = 100*mean(tn_por_exame ,'all')
% 
% tp_por_exame = squeeze(any(tester.TP(tester.dataloader.signalFrequencies,:,:)>0,2));
% tp_medio = 100*mean(tp_por_exame ,'all')
% 
% fn_por_exame = squeeze(any(tester.FN(tester.dataloader.signalFrequencies,:,:)>0,2));
% fn_medio = 100*mean(fn_por_exame ,'all')

% find non-empty parameter sets


fp_por_param = cell(1);
tp_por_param = cell(1);
fn_por_param = cell(1);
tn_por_param = cell(1);

for params_idx = 1:numel(selected_exams)
for stim = tester.stimulusIndices
    for subj = tester.subjectIndices
        selected_exams = tester.groupDecisions(1, 1,:);

        
            if ~isempty(cell2mat(selected_exams(params_idx)))

                decision_group_idx = sub2ind(size(tester.groupDecisions), ...
                    stimulusIndex, subjectIndex, params_idx );

                p.nWindows = cell2mat(obj.nWindows(decision_group_idx));
                p.K_stages = cell2mat(obj.K_stages(decision_group_idx));
                
                M = p.nWindows;
                K = p.K_stages;
                
                obj = obj.single_exam_beta_cgst_threshold(M,K);
                
                obj.groupStageAlphas{epoch_idx} = obj.stageAlphas;
                obj.groupStageGammas{epoch_idx} = obj.stageGammas;

            end
        end
    end
end


some_stim = ordc.stimulusIndices(end);
some_subj = 1; % 1 - 11
this_msc = cell2mat(ordc.groupMSC(5,5,K_stages));
this_msc = this_msc(:,:,single_channel);
for k=1:size(this_msc,2)
    s1 = stem(sum(this_msc(:,1:k,1),2));
    hold on
    title(['Epoch = ',num2str(k)])
    s2 = stem(tester.signalFrequencies, ...
        sum(this_msc(tester.signalFrequencies,1:k,1),2), 'filled', 'color', 'k');

    s3 = stem(tester.noiseFrequencies, ...
        sum(this_msc(tester.noiseFrequencies,1:k,1),2), 'filled', 'color', 'r');

    l1 = yline(tester.stageAlphas(k), 'g--');
    l2 = yline(tester.stageGammas(k), 'r--');

    legend('MSC', 'Signal', 'Noise', 'Detect above this', 'Stop below this')
    grid on
    xlim([0,501])
    hold off
    pause(0.66*2.2)
end

% 
% % Show 3 exam timeseries and PSD (from random random stimulus-subject-electrode)
% randomly_shown_subj = dtl.selectedZanoteliSubjects(randperm(nSubj,3));
% randomly_shown_stim = dtl.selectedZanoteliStimuli(randi([1 2],1,3));
% randomly_shown_electrode = [single_channel single_channel single_channel]; %randperm(dtl.nChannels,3)
% 
% % find nonempty epoch parameters
% nonempty_params = [];
% for i = 1:numel(ordc.epochs(dtl.selectedZanoteliStimuli(1),dtl.selectedZanoteliSubjects(1),:))
%     if ~isempty(cell2mat(ordc.epochs(dtl.selectedZanoteliStimuli(1),dtl.selectedZanoteliSubjects(1),i)))
%         nonempty_params = [i nonempty_params];
%     end
% end
% randomly_selected_epoch_param = nonempty_params(randi([1 numel(nonempty_params)],1,3));
% 
% 
% figure(1)
% for i = 1:3
% 
%     current_signal = cell2mat( ...
%         dtl.groupSignals(randomly_shown_stim(i),randomly_shown_subj(i)));
% 
%     current_SIGNAL =  cell2mat( ...
%         ordc.groupMSC(randomly_shown_stim(i),randomly_shown_subj(i),randomly_selected_epoch_param(i)));
% 
%     random_epoch = randi([1,size(current_signal,2)],1);
%     random_test = randi([1,size(current_SIGNAL,2)-1],1);
% 
%     current_signal = current_signal(:,random_epoch,randomly_shown_electrode(i));    
%     lead_name = dtl.zanoteliLeads(randomly_shown_electrode(i));
% 
%     % add 2 seconds removed during preprocessing
%     exam_time = ((random_epoch+2)*dtl.fs:(random_epoch+2)*dtl.fs+numel(current_signal)-1)';
% 
%     subplot(2,3,i)
%     plot(exam_time,current_signal)
%     grid on 
%     xlabel('Time [ms]')
%     ylabel('Voltage [V]')
%     title(['Measures Aquired on' lead_name])
% 
%     subplot(2,3,i+3)
%     stem(abs(current_SIGNAL(:,random_test,randomly_shown_electrode(i))), ...
%         'filled', 'MarkerSize',3,'LineWidth',0.1)
%     hold on
%     stem(dtl.signalFrequencies, ...
%         abs(current_SIGNAL(dtl.signalFrequencies,random_test,randomly_shown_electrode(i))), ...
%         'filled', 'MarkerSize',3,'LineWidth',0.1, 'Color','red');
%     grid on 
%     xlabel('Frequency [Hz]')
%     ylabel('MSC')
%     title(['MSC computed from' lead_name])
%     hold off
% 
% end
