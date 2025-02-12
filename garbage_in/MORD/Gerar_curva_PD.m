
clear all
close all
clc


det = {'MMSC';'aMSC';'pMSC';'MCSM';'aCSM';'pCSM';'MLFT';'aLFT';'pLFT';'MGBT';'aGBT';'pGBT'};

% methods = ['time_Cholesky_corrected'];
methods = ['Monte_Carlo_default'];
Mmax = 10;
%A = [1 .3;.7 1];
A = [1 0;0 1];

for ii = 1:9    
       detector = cell2mat(det(ii,:));
       PDtime(detector,methods,Mmax,A) 
end

%% 
hPD =[];
for ii = 1:9 
       detector = cell2mat(det(ii,:));
%       load([ detector '_' methods '_' num2str(Mmax) '_filter'],'PD','vetor_SNR','A','parameters');
        load([ detector '_' methods '_' num2str(Mmax)],'PD','vetor_SNR','A','parameters');    
        hPD(:,ii) = PD 
end


