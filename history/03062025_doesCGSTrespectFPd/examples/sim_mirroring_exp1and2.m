%% Clear workspace and initialize parameters
clearvars; close all; clc;

%% Test pipeline
% Apply TEST(s) to objective response detector (MSC)
% tester = ORDTester(ordc);
% tester = tester.compute_beta_cgst_thresholds();

tester = ORDTester(ORDCalculator(DataLoader('sim')));

% tester.dataloader = tester.dataloader.resetDuration(55);
% tester.dataloader = tester.dataloader.genBulkSim( ...
%     groupNoiseMean=[-30 -20], groupNoiseStd=5*randi(5,1,11)).computeBulkFFTs(mode='sim');
% tester.dataloader = tester.dataloader.resetDuration(90);
% tester.dataloader = tester.dataloader.genBulkSim( ...
%     groupNoiseMean=[-20], groupNoiseStd=randi(5,1,11)).computeBulkFFTs(mode='sim');
% tester.age() %                                                                                      
tester.dataloader = tester.dataloader.resetDuration(90);
tester.dataloader = tester.dataloader.genBulkSim( ...
    groupNoiseMean=[-30 -28 -32], groupNoiseStd=randi(5,1,11)).computeBulkFFTs(mode='sim');
tester.age() %   
% ppc = PreProcessor();
% ppc = ppc.bulkZanoteliPreprocess(tester.dataloader);
% ppc = ppc.bulkAntunesFilter(tester.dataloader);
% dtl.groupSignals = ppc.groupProcessedSignals;
% tester.dataloader.groupSignals = ppc.groupFilteredSignals;

% tester.dataloader = tester.dataloader.computeBulkFFTs(mode='sim');

%%
% Pretty much instantly runs
% K_stages = 5;
% tester.ord_calculator = ORDCalculator(tester.dataloader).fit_epochs( ...
%     stimulusIndices = tester.dataloader.groupNoiseMean,...
%     subjectIndices = tester.dataloader.groupNoiseStd,...
%     startWindows = 1, ...
%     K_stages = K_stages, ...
%     single_or_bulk = 'bulk',...
%     lastWindowCalcMethod = 'exactK', ... % maxFromStart, maxFromLast, exactK, fromSizeType
%     sizeType = 'fixedSize'... % minToMax, minToFix, withResampling, default = fixedSize
%     ... % then, compute on selected channels:
%     );

% Takes at least about 860.52 seconds = 15 mins to run 
% K_stages = [2 3 5 8 10];
tester.ord_calculator = ORDCalculator(tester.dataloader).fit_epochs( ...
    stimulusIndices = tester.dataloader.groupNoiseMean,...
    subjectIndices = tester.dataloader.groupNoiseStd,...
    startWindows = [1:5 20 30], ... % windowStepSizes = 52, ...
    K_stages = 2:10,...
    single_or_bulk = 'bulk',...
    lastWindowCalcMethod = 'exactK', ... % maxFromStart, maxFromLast, exactK, fromSizeType
    sizeType = 'fixedSize'... % minToMax, minToFix, withResampling, default = fixedSize
    ... % then, compute on selected channels:
    );

fprintf('\t [%s] Epochs are computed.\n', datetime())
tester.ord_calculator.age()

%%

single_channel = 1:16;
tester.ord_calculator = tester.ord_calculator.bulk_compute_msc(channels = single_channel);
tester.age()

%% Test pipeline
% Apply TEST(s) to objective response detector (MSC)
tester = ORDTester(tester.ord_calculator);
tester.desired_alpha = 0.05;
tester = tester.compute_bulk_sim_beta_cgst_decisions();
tester.age()

%% Show results


nParams = numel(tester.epochs(1,1,:));
nChannels = 16; %numel(tester.ord_calculator.channels);
tester.FP = zeros(nParams,nChannels);
tester.FN = zeros(nParams,nChannels);
tester.TP = zeros(nParams,nChannels);
tester.TN = zeros(nParams,nChannels);

stimulusIndices = 1:numel(tester.ord_calculator.dataloader.groupNoiseMean);
subjectIndices = 1:numel(tester.ord_calculator.dataloader.groupNoiseStd);

for stimulusIndex = stimulusIndices
    for subjectIndex = subjectIndices
        for channel_idx = 1:numel(tester.ord_calculator.channels)
            nonemptyparams_idxs = [];
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
                    
                    test_channel_idx =tester.ord_calculator.channels(channel_idx);
                    
                    tester.FN(epoch_param_idx, test_channel_idx) = ...
                            sum(any(exam_fn(tester.signalFrequencies,:,test_channel_idx)>0,2))...
                            + tester.FN(epoch_param_idx, test_channel_idx);

                    tester.TP(epoch_param_idx, test_channel_idx) = ...
                            sum(any(exam_tp(tester.signalFrequencies,:,test_channel_idx)>0,2))...
                            + tester.TP(epoch_param_idx, test_channel_idx);

                    tester.FP(epoch_param_idx, test_channel_idx) = ...
                            sum(any(exam_fp(tester.noiseFrequencies,:,test_channel_idx)>0,2))...
                            + tester.FP(epoch_param_idx, test_channel_idx);

                    tester.TN(epoch_param_idx, test_channel_idx) = ...
                            sum(any(exam_tn(tester.noiseFrequencies,:,test_channel_idx)>0,2))...
                            + tester.TN(epoch_param_idx, test_channel_idx);

                    nonemptyparams_idxs = [nonemptyparams_idxs,...
                                            epoch_param_idx];
    
                end
            end
        end
    end
end


denom = numel(tester.noiseFrequencies)*numel(tester.ord_calculator.dataloader.groupNoiseMean)...
    *numel(tester.ord_calculator.dataloader.groupNoiseStd);
fp_rate = tester.FP/(denom);

tp_rate = tester.TP/(denom);

fn_rate = tester.FN/(denom);

tn_rate = tester.TN/(denom);

confmat = table( ...
    100*mean(fn_rate(nonemptyparams_idxs,single_channel),'all'), ...
    100*mean(fp_rate(nonemptyparams_idxs,single_channel),'all'), ...
    100*mean(tp_rate(nonemptyparams_idxs,single_channel),'all'), ...
    100*mean(tn_rate(nonemptyparams_idxs,single_channel),'all'), ...
    'VariableNames',{'fn','fp','tp','tn'})
% 
% tester_checkpoint = struct();
% tester_checkpoint.epochs = tester.epochs ;
% tester_checkpoint.groupTP = tester.groupTP;
% tester_checkpoint.groupTN = tester.groupTN;
% tester_checkpoint.groupFP = tester.groupFP;
% tester_checkpoint.groupFN = tester.groupFN;
% tester_checkpoint.previous_cgst_thresholds = tester.previous_cgst_thresholds;
% tester_checkpoint.allTestFrequencies = tester.allTestFrequencies;
% 
% tester_checkpoint.selectedZanoteliSubjects = tester.dataloader.selectedZanoteliSubjects;
% tester_checkpoint.selectedZanoteliStimuli = tester.dataloader.selectedZanoteliStimuli;
% tester_checkpoint.epochs_index_metadata = tester.ord_calculator.epochs_index_metadata;
% 
% tester_checkpoint.vaz_translated_Kstages = vaz_translated_Kstages;
% tester_checkpoint.vaz_startWindows = vaz_startWindows;
% 

[Y,MO,D,H,MI,S] = datevec(datetime);
s = ['tester_checkpoint_1_T',...
    num2str(H),'_',num2str(MI),'_',num2str(fix(S)),...
    num2str(D),'_',num2str(MO),'_',num2str(Y),'.mat']
% save(s,'tester_checkpoint','-v7.3')


%%
figure(1)
subplot(211)
plot(100*fp_rate(nonemptyparams_idxs,single_channel),'.')
subplot(212)
plot(100*tp_rate(nonemptyparams_idxs,single_channel),'.')

figure(2)
subplot(211)
b1 = boxchart(100*fp_rate(nonemptyparams_idxs,single_channel));
subplot(212)
b2 = boxchart(100*tp_rate(nonemptyparams_idxs,single_channel));

% save(s,'-v7.3')
