%% Clear workspace and initialize parameters
clearvars; close all; clc;

%% Test pipeline
% Apply TEST(s) to objective response detector (MSC)
% tester = ORDTester(ordc);
% tester = tester.compute_beta_cgst_thresholds();

tester = ORDTester(ORDCalculator(DataLoader('sim')));
tester.desired_alpha = 0.05;

snr = -30
[params,tester ]= tester.validateDetectionThresholds(SNRmean = snr, duration = 90);


check = tester.FP+tester.FN+tester.TP+tester.TN;
check_tested = check(params.allTestFrequencies,:,:);
check_untested = check;
check_untested(params.allTestFrequencies,:,:) = [];

% decisao_por_exame = zeros(params.dataloader.noiseFrequencies,size(tester.FP,3));
fp_por_exame = squeeze(any(tester.FP(params.dataloader.noiseFrequencies,:,:)>0,2));
fp_medio = 100*mean(fp_por_exame ,'all') 

tn_por_exame = squeeze(any(tester.FP(params.dataloader.noiseFrequencies,:,:)>0,2));
tn_medio = 100*mean(tn_por_exame ,'all') 

tp_por_exame = squeeze(any(tester.TP(params.dataloader.signalFrequencies,:,:)>0,2));
tp_medio = 100*mean(tp_por_exame ,'all') 

fn_por_exame = squeeze(any(tester.FN(params.dataloader.signalFrequencies,:,:)>0,2));
fn_medio = 100*mean(fn_por_exame ,'all') 



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
                exam_fp = cell2mat(tester.groupFP(stimulusIndex, subjectIndex, epoch_param_idx));
                if ~isempty(exam_fp)
                    % If there was a detection on any stage, 
                    % assign 1 to that exam, and repeat for all
                    % frequencies. Then, sum all 1s.
                    % + Add previous results from other subjs
                    % To compute rate (pct), divide by #freqs and #subj

                    exam_fn = cell2mat(tester.groupFN(stimulusIndex, subjectIndex, epoch_param_idx));
                    exam_tp = cell2mat(tester.groupTP(stimulusIndex, subjectIndex, epoch_param_idx));
                    exam_tn = cell2mat(tester.groupTN(stimulusIndex, subjectIndex, epoch_param_idx));

                    tester.FP(epoch_param_idx, channel_idx) = ...
                            sum(any(exam_fp(tester.noiseFrequencies,:,channel_idx)>0,2))...
                            + tester.FP(epoch_param_idx, channel_idx);

                    tester.FN(epoch_param_idx, channel_idx) = ...
                            sum(any(exam_fp(tester.signalFrequencies,:,channel_idx)>0,2))...
                            + tester.FP(epoch_param_idx, channel_idx);

                    tester.TP(epoch_param_idx, channel_idx) = ...
                            sum(any(exam_fp(tester.signalFrequencies,:,channel_idx)>0,2))...
                            + tester.FP(epoch_param_idx, channel_idx);

                    tester.TN(epoch_param_idx, channel_idx) = ...
                            sum(any(exam_fp(tester.noiseFrequencies,:,channel_idx)>0,2))...
                            + tester.FP(epoch_param_idx, channel_idx);
    
                end
            end
        end
    end
end


fp_rate = tester.FP/(numel(tester.noiseFrequencies)+numel(tester.subjectIndices))

tp_rate = tester.TP/(numel(tester.noiseFrequencies)+numel(tester.subjectIndices))

fn_rate = tester.FN/(numel(tester.noiseFrequencies)+numel(tester.subjectIndices))

tn_rate = tester.TN/(numel(tester.noiseFrequencies)+numel(tester.subjectIndices))
