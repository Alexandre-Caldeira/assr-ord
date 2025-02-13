%%
 load(['timeM_50db_M1' '.mat'],'timeM')
load(['pareto_50db_M1' '.mat'],'TXD','Mmax','parametros');
%%
load(['timeM_50db_M4' '.mat'],'timeM')
load(['pareto_50db_M4' '.mat'],'TXD','Mmax','parametros');

%%
% taxa de detecção x tempo 
% timeM = time(binsM,:); 
% timeM(timeM==-1) = Mmax;
% timeM = mean(timeM,1)'*1; %1segundo por janela
% TXD = mean(mean(Tdr(binsM,:,:),3),1)' *100;


% figure1 = figure;
axes1 = axes('Parent',figure1); hold(axes1,'on');  
plot([0 1]*100, [Mmax Mmax],'-.k','linewidth',1) 
plot([TXD(end) TXD(end)], [min(timeM) max(timeM)],'-.k','linewidth',1) 


for ii = 1:size(parametros,1)
    plot(TXD(ii),timeM(ii),'.k','Markersize',6,'DisplayName',[num2str(parametros(ii,1)) '-' num2str(parametros(ii,2))])
end

[ p, idxs] = paretoFront([TXD,(-timeM)] ); 
auxL = p(:,1)<0.5; 
p(auxL,:) = [];
idxs(auxL,:) = [];

[~,ind] = sort(p(:,1));
p = p(ind,:);
idxs = idxs(ind,:);
% p([3,4,7,9],:) = [];
% idxs([3,4,7,9],:) = [];

plot(TXD(idxs),timeM(idxs),'-ok','Markersize',8,'linewidth',1.2) 

% ylim([min(-p(:,2))*.9 max(-p(:,2))*1.1])
% xlim([min(p(:,1))*80 max(p(:,1))*104])

set(axes1,'XMinorTick','on');
set(axes1,'YMinorTick','on');
box(axes1,'off');




for ii = 1:size(idxs,1)

    [I] = find((TXD == TXD(idxs(ii))) & (timeM==timeM(idxs(ii))));
    
    fprintf('\nPD = %f Tempo = %f ',TXD(idxs(ii)),timeM(idxs(ii))); 
    fprintf(' NI = %d ',length(I)); 
     I = I(1); 
    
    for jj = I  
        
        fprintf(' - Buffer:%d, M_step:%d', parametros(jj,1),parametros(jj,2)); 
    %    text(PD_T(jj)*99,T_exame_G(idxs(ii))*.985,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) ','  num2str(parametros(jj,3)) '\}' ]);
%            text(PD_T(jj)*99,T_exame_G(idxs(ii))*.985,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) ',100,'  num2str(parametros(jj,3)) '\}' ]);
%        text(PD_T(jj)*99,T_exame_G(idxs(ii))*.985,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) ',100\}' ]);
       text(TXD(jj),timeM(idxs(ii))*.975,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) '\}' ]);
    end    
end
fprintf('\n'); 
xlim([min(TXD(idxs))*.95,max(TXD(idxs))*1.05])
ylim([min(timeM(idxs))*.95 Mmax*1.05])

grid on
hold on

% %%
% hold on
% %%
% xlabel('Taxa de Detecção (%)','fontsize',12); 
% ylabel('Tempo Médio de Exame (s)','fontsize',12);
% title(['Estímulo ', (Intensidade{1}), ': Método 04'])
% % title(['Estímulo ', (Intensidade{1}), ' SPL: Método 04'])

