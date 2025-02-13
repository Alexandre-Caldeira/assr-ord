%function analisar_resultadosv3()
% clear all
%close all
clc 

mord = {'MMSC';'aMSC';'pMSC';'MCSM';'aCSM';'pCSM';'MLFT';'aLFT';'pLFT'};
% mord = {'aMSC'; 'aCSM';'aLFT';'pMSC';'pCSM';'pLFT';'MMSC';'MCSM';'MLFT'}
% mord = { 'MMSC'; 'aMSC'; 'pMSC';'MCSM';'aCSM';'pCSM';'MLFT';'aLFT';'pLFT';'MGBT';'aGBT';'pGBT'};
% detS= { '$\hat{\kappa}_{N}^{2}(f)$'; '$\overline{\hat{\kappa}_{y}^{2}}(f)$'; '$\sqrt[N]{\hat{\kappa}_{N!}^{2}(f)}$';...}
%     '$\hat{\rho}_{N}^{2}(f)$'; '$\overline{\hat{\rho}_{y}^{2}}(f)$'; '$\sqrt[N]{\hat{\rho}_{N!}^{2}(f)}$';...
%     '$\hat{\phi}_{N}(f_{o})$'; '$\overline{\hat{\phi}_{y}}(f_{o})$'; '$\sqrt[N]{\hat{\phi}_{N!}(f_{o})}$';...
%  '$\hat{\beta}_{N}(f)$'; '$\overline{\hat{\beta}_{y}}(f)$'; '$\sqrt[N]{\hat{\beta}_{N!}(f)}$'};
detS= { '$\hat{\kappa}_{N}^{2}(f)$'; '$\overline{\hat{\kappa}_{y}^{2}}(f)$'; '$\sqrt[N]{\hat{\kappa}_{y}^{2}(f)}$';...}
    '$\hat{\rho}_{N}^{2}(f)$'; '$\overline{\hat{\rho}_{y}^{2}}(f)$'; '$\sqrt[N]{\hat{\rho}_{y}^{2}(f)}$';...
    '$\hat{\phi}_{N}(f_{o})$'; '$\overline{\hat{\phi}_{y}}(f_{o})$'; '$\sqrt[N]{\hat{\phi}_{y}(f_{o})}$';...
 '$\hat{\beta}_{N}(f)$'; '$\overline{\hat{\beta}_{y}}(f)$'; '$\sqrt[N]{\hat{\beta}_{y}(f)}$'};


 %Limite de confiança --------------------------
nh = 20; %número de harmônicos no Fp 
nv = 10; %número de voluntários 
NN = nh *nv; 
p = 0.05; %nível de significãncias 
limite_conf_sup = .90; %intervalo de confiança
%se aumentar para 0.95 -> para de 90% limite de 92%
%se aumentar para 0.90 -> limite em simulação de 85%.
pp = (1-limite_conf_sup)/2; 

Nrun = 100000;
for ii = 1:Nrun 
r(ii) = binornd(NN,p);
end

ni = quantile(r,1-limite_conf_sup);
ns = quantile(r,limite_conf_sup);
mean(r>=ni);
mean(r<=ns);
mean((r<=ns).*(r>=ni));


% linf = (quantile(r,1-limite_conf_sup))/NN;
% lsup = (quantile(r,limite_conf_sup))/NN;
% 
% lsup + linf;

linf = (quantile(r,pp)+1)/NN;
lsup = (quantile(r,1-pp)-1)/NN;

n_sinais = [1 2 4 6 8];


%% 
figure
lw =1.2;
cont = 1; 
for ii = 1:9

    det = cell2mat(mord(ii));
    load([det '_vectors.mat'])  
    
    subplot(3,3,cont) 
    
 % --------------------------------------------------------   
%     Nrun = 100000;
%     p = eval(['DRMC_' det '(1);'])/100;
%     for kk = 1:Nrun 
%     r(kk) = binornd(NN,p);
%     end
% 
%     ni = quantile(r,1-limite_conf_sup);
%     ns = quantile(r,limite_conf_sup);
%     mean(r>=ni);
%     mean(r<=ns);
%     mean((r<=ns).*(r>=ni));
% 
% 
%     linf = (quantile(r,1-limite_conf_sup))/NN;
%     lsup = (quantile(r,limite_conf_sup))/NN;    
% --------------------------------------------------  
%     rectangle('Position',... 
%           [1 linf*100 8 (lsup-linf)*100],...
%           'FaceColor',[0.5 1 0.5 0.5],'EdgeColor','none')
 %area([0 9], [lsup*100 lsup*100],linf*100,'FaceColor',[.8 .8 .8])

    hold on 
    errorbar(n_sinais,eval(['DRMC_' det]),eval(['stdDRMC_' det]),'-k','LineWidth',lw);
    hold on
        if ((ii)>=7)
    xlabel('Number of signals') ;
    end
    
    if (((ii)==4))
        ylabel(['Detection/ False Positive rates (%)'], 'fontsize',14)
    end 
    
    
    
    

 %   errorbar(n_sinais,eval(['DRMC_CORR_' det]),eval(['stdDRMC_CORR_' det]),'--sk','LineWidth',lw);
 %   plot(n_sinais,eval(['DRMC_CORR_' det]),'-*k','LineWidth',lw);
 
%------------------------------------------------
    r=[];
    Nrun = 100000;
    p = eval(['FPMC_' det '(1);'])/100;
    for kk = 1:Nrun 
    r(kk) = binornd(NN,p);
    end
    pp = (1-limite_conf_sup)/2; 
    
    ni = (quantile(r,pp));
    ns = (quantile(r,1-pp));
    mean(r>=ni);
    mean(r<=ns);
    mean((r<=ns).*(r>=ni))

%    linf = (quantile(r,pp)+1)/NN;
     linf = (quantile(r,pp))/NN;
 %   lsup = (quantile(r,1-pp)-1)/NN;  
    lsup = (quantile(r,1-pp))/NN; 
%--------------------------------------------------  
area([0 9], [lsup*100 lsup*100],linf*100,'FaceColor',[.8 .8 .8])
 
load([det '_CORR_vectors_VCESP_MediaCoerente.mat'])

errorbar(n_sinais,eval(['FPMC_' det]),eval(['stdFPMC_' det]),'--k','LineWidth',lw);

 errorbar(n_sinais,eval(['DRMC_CORR_' det]),eval(['stdDRMC_CORR_' det]),'b','LineWidth',lw);
    errorbar(n_sinais,eval(['FPMC_CORR_' det]),eval(['stdFPMC_CORR_' det]),'--b','LineWidth',lw); 
 
    if (ii<10)
        %load([det '_CORR_vectors_VCESP_MediaCoerente.mat'])
        %errorbar(n_sinais,eval(['DRMC_CORR_' det]),eval(['stdDRMC_CORR_' det]),'r','LineWidth',lw);
        %errorbar(n_sinais,eval(['FPMC_CORR_' det]),eval(['stdFPMC_CORR_' det]),'--r','LineWidth',lw); 
    
        load([det '_CORR_vectors_freq_VCEST.mat'])
        errorbar(n_sinais,eval(['DRMC_CORR_' det]),eval(['stdDRMC_CORR_' det]),'r','LineWidth',lw);
        errorbar(n_sinais,eval(['FPMC_CORR_' det]),eval(['stdFPMC_CORR_' det]),'--r','LineWidth',lw); 
    
    end
    
%     load([det '_CORR_vectors_freq_VCEST.mat'])
%     errorbar(n_sinais,eval(['DRMC_CORR_' det]),eval(['stdDRMC_CORR_' det]),'m','LineWidth',lw);
%      errorbar(n_sinais,eval(['FPMC_CORR_' det]),eval(['stdFPMC_CORR_' det]),'m','LineWidth',lw); 
    
    
    hold off
 %   legend({'Using defaut critical values','Using Cholesky-corrected critical values'},'edgecolor','none')
    axis([1 8.1 0 80])
    %axis([1 8.1 0 100])
      title([detS(ii)],'Interpreter','latex','fontsize',12)

      
    cont = cont+1;
end



%% plotar a diferença entre os dois valores críticos 
% figure
% lw =1.2;
% cont = 1; 
% for ii = 1:9
% 
%     det = cell2mat(mord(ii));
%     load([det '_CORR_vectors_VCESP.mat'])
%     subplot(3,3,cont) 
%     DR1 = eval(['DRMC_CORR_' det]);
%     FP1 = eval(['FPMC_CORR_' det]);
%     load([det '_CORR_vectors_VCEST.mat'])
%     DR2 = eval(['DRMC_CORR_' det]);
%     FP2 = eval(['FPMC_CORR_' det]);
%     
%     plot(n_sinais,(DR2-DR1),'k','LineWidth',lw)
%     hold on 
%     plot(n_sinais,(FP2-FP1),'r','LineWidth',lw)
%     plot([0 10],[0 0],'-b','LineWidth',lw)
%     
%     hold off 
%     
% 
%     if ((ii)>=7)
%      xlabel('Number of signals') ;
%     end
%     
%     if (((ii)==4))
%         ylabel(['Detection/ False Positive rates (%)'], 'fontsize',14)
%     end 
%     
%     
%    axis([2 8.1 -5 5])
%    title([detS(ii)],'Interpreter','latex','fontsize',12)  
%     cont = cont+1;
% end


%%  Maiores taxas de detecção 
n = [1 2 4 6 8];
TXDm = []; 
TXD_sem =[];
FP_sem=[];
Nbest =[]; 
limiar = {'maior';'menor';'igual'} ;
for ii = 1:9
    det = cell2mat(mord(ii));
    load([det '_CORR_vectors_VCEST_MediaCoerente.mat']);
   % load([det '_vectors_PCA.mat'])
    TXD = eval(['DRMC_CORR_' det]);
    FP = eval(['FPMC_CORR_' det]);
  
    [~,a] = max(TXD);
    TXDm(ii,1) = TXD(a);
    Nbest(ii,1) = n(a);
    FPbest(ii,1) = FP(a);
    %txv(1:10,ii) = DRMC_CORR(:,a);
    
    if FPbest(ii,1)>(lsup*100) 
        inf(ii,1) = limiar(1);
    elseif FPbest(ii,1)<(linf*100)
        inf(ii,1) = limiar(2);
    else
        inf(ii,1) = limiar(3);
    end
    
     load([det '_vectors.mat']);
    %equivalente sem correção 
    TXD = eval(['DRMC_' det]);
    FP = eval(['FPMC_' det]);
    TXD_sem(ii,1) = TXD(a);
    FP_sem(ii,1) = FP(a);    
end

[~,a] = sort(TXDm,'descend');
tt=[];
tt = table(cell2mat(mord(a)),TXDm(a,1),Nbest(a,1),FPbest(a,1),cell2mat(inf(a)),TXD_sem(a,1),FP_sem(a,1)) 

txv = txv(:,a');

%compare 
[p,t,stats] = anova1(txv);

[c,m,h,nms] = multcompare(stats);

[~,a] = sort(TXDm);
% Set the remaining axes properties
set(gca,'YTick',1:12,'YTickLabel',...
    mord(a));

% %% tempo 
% lw =1.5;
% nRuns = 70000;
% for ii = 1:12
%     det = cell2mat(mord(ii));
%     load([det '_vectors.mat']);
%     tempoMC(ii,1) = eval(['tempoMC_' det])/nRuns;
%     tempoChol(ii,1) = eval(['tempoChol_' det])/nRuns;
% end
% 
% tempoChoMedio = (tempoChol- tempoMC)/(nv*2);
% %tempoChoMedio = (tempoChol- tempoMC)/(1);


% [~,a] = sort(tempoChoMedio);
% figure
% plot(1:12, tempoMC(a),'k','LineWidth',lw);ylabel('Time (s)')
% hold on
% plot(1:12, tempoChoMedio(a),'--k','LineWidth',lw);ylabel('Time (s)')
% legend({'Using defaut critical values','Using Cholesky-corrected critical values'},'edgecolor','none','fontsize',10)

% set(gca,'XTick',1:12,'XTickLabel',...
%     mord(a),'fontsize',9);
% set(gca,'XTick',1:12,'XTickLabel',...
%     detS2(a),'fontsize',8);
% set(gca,'XTick',1:12,'XTickLabel',...
%     detS3(a),'fontsize',12);
% 
% set(gca,'TickLabelInterpreter','latex')
% xlim([1 12])




%% aumento no FP 
% for ii = 1:12
% 
%     det = cell2mat(mord(ii));
%     load([det '_vectors.mat'])
%     x1 = eval(['fpMC_' det]);
%     x2 = eval(['FPMC_CORR_' det]);
%     
%     y1(:,ii) = x1./(x1(1));
%     y2(:,ii) = x2./(x2(1));
% end
% 
% figure
% plot(y1,'-')
% hold on 
% plot(y2,':')

%%






% for ii = 1:size(mord,1)
% 
%     det = cell2mat(mord(ii));
%     load([det '_vectors.mat'])
% 
%    
%     figure;
%     
%     yyaxis left
%     errorbar(n_sinais,eval(['DRMC_' det]),eval(['stdDRMC_' det]),'b');xlabel('Number of signals');ylabel('Detection rate (%)')
%     hold on
%     errorbar(n_sinais,eval(['DRMC_CORR_' det]),eval(['stdDRMC_CORR_' det]),'--b');xlabel('Number of signals');ylabel('Detection rate (%)')
%     axis([.7 8.3 0 80])
%     yyaxis right
%     area([0 9], [lsup*100 lsup*100],linf*100,'LineStyle',':','FaceColor',[.8 .8 .8])
%     errorbar(n_sinais,eval(['fpMC_' det]),eval(['stdfpMC_' det]),'r');xlabel('Number of signals');ylabel('False Detection rate (%)')
%     errorbar(n_sinais,eval(['FPMC_CORR_' det]),eval(['stdFPMC_CORR_' det]),'--r');xlabel('Number of signals');ylabel('False Detection rate (%)')
% %     plot([0 9],[5 5],':k')
% %     plot(0 9,  [linf linf]., 
%     hold off
%     legend({'Using defaut critical values','Using Cholesky-corrected critical values'},'edgecolor','none')
%     axis([.7 8.3 0 80])
%     title(det)
% 
% end
%     
%end

    
    
    
