
figure1 = figure;
axes1 = axes('Parent',figure1); 
hold(axes1,'on');  

% Single Shot %
%plot([TXD_2006(1) TXD_2006(1)], [min(timeM_2006) max(timeM_2006)],'-.k','linewidth',1) 
%ERRO
%plot([TXD_2013(1) TXD_2013(1)], [min(timeM_2013) max(timeM_2013)],'-.k','linewidth',1) 
%ERRO
%plot([TXD_2005(1) TXD_2005(1)], [min(timeM_2005) max(timeM_2005)],'-.k','linewidth',1) 
%ERRO

% PLOTAR OS PONTOS %
%for ii = 1:(size(parametros_par_2006))
%    plot(TXD_2006(ii),timeM_2006(ii),'.k','Markersize',6,'DisplayName',[num2str(parametros_par_2006(ii,1)) '-' num2str(parametros_par_2006(ii,2))])
%end
%ERRO
%for ii = 1:(size(parametros_par_2013))
%    plot(TXD_2013(ii),timeM_2013(ii),'.k','Markersize',6,'DisplayName',[num2str(parametros_par_2013(ii,1)) '-' num2str(parametros_par_2013(ii,2))])
%end
%ERRO
%for ii = 1:(size(parametros_par_2005))
%    plot(TXD_2005(ii),timeM_2005(ii),'.k','Markersize',6,'DisplayName',[num2str(parametros_par_2005(ii,1)) '-' num2str(parametros_par_2005(ii,2))])
%end

%ERRO

%TXD_2005 = TXD;
%timeM_2005 = timeM;
%FP_2005 = FP;
%parametros_par_2005 = parametros_par;


[ p, idxs] = paretoFront([TXD_2005,(-timeM_2005)] ); 
auxL = p(:,1)<0.5; 
p(auxL,:) = [];
idxs(auxL,:) = [];
[~,ind] = sort(p(:,1));
p = p(ind,:);
idxs = idxs(ind,:);
idxs_2005 = idxs;
% p([3,4,7,9],:) = [];
% idxs([3,4,7,9],:) = [];
plot(TXD_2005(idxs_2005),timeM_2005(idxs_2005),'-ob','Markersize',8,'linewidth',1.8) 
%ERRO



%ERRO

%TXD_2006 = TXD;
%timeM_2006 = timeM;
%FP_2006 = FP;
%parametros_par_2006 = parametros_par;

[ p, idxs] = paretoFront([TXD_2006,(-timeM_2006)] ); 
auxL = p(:,1)<0.5; 
p(auxL,:) = [];
idxs(auxL,:) = [];
[~,ind] = sort(p(:,1));
p = p(ind,:);
idxs = idxs(ind,:);
% p([3,4,7,9],:) = [];
% idxs([3,4,7,9],:) = [];
idxs_2006 = idxs;
plot(TXD_2006(idxs_2006),timeM_2006(idxs_2006),'-oc','Markersize',8,'linewidth',1.8) 


%ERRO

%TXD_2013 = TXD;
%timeM_2013 = timeM;
%parametros_par_2013 = parametros_par;

%[ p, idxs] = paretoFront([TXD_2013,(-timeM_2013)] ); 
%auxL = p(:,1)<0.5; 
%p(auxL,:) = [];
%idxs(auxL,:) = [];
%[~,ind] = sort(p(:,1));
%p = p(ind,:);
%idxs = idxs(ind,:);
%idxs_2013 = idxs;
% p([3,4,7,9],:) = [];
%plot(TXD_2013(idxs_2013),timeM_2013(idxs_2013),'-ob','Markersize',8,'linewidth',1.2) 


%ERRO

%TXD_2015 = TXD;
%timeM_2015 = timeM;
%FP_2015 = FP;
%parametros_par_2015 = parametros_par;


[ p, idxs] = paretoFront([TXD_2015,(-timeM_2015)] ); 
auxL = p(:,1)<0.5; 
p(auxL,:) = [];
idxs(auxL,:) = [];
[~,ind] = sort(p(:,1));
p = p(ind,:);
idxs = idxs(ind,:);
idxs_2015 = idxs;
% p([3,4,7,9],:) = [];
% idxs([3,4,7,9],:) = [];
plot(TXD_2015(idxs_2015),timeM_2015(idxs_2015),'-om','Markersize',8,'linewidth',1.8) 
%ERRO

set(axes1,'XMinorTick','on');
set(axes1,'YMinorTick','on');
box(axes1,'off');

axes1.XColor = 'k';  % Define a cor do eixo X
axes1.YColor = 'k';  % Define a cor do eixo Y

% Configuração das linhas pontilhadas para os eixos X e Y
%axes1.XAxis.LineStyle = '--';  % Define o estilo da linha do eixo X como pontilhada
%axes1.YAxis.LineStyle = '--';  % Define o estilo da linha do eixo Y como pontilhada

hold off
xlabel('Detection Rate (%)','fontsize',15); 
ylabel('Detection Time(s)','fontsize',15);
fprintf('\n'); 


hold on
plot([0 1]*100, [240 240],'-.k','linewidth',1) 


fprintf('\n'); 
xlim([53, 69])
ylim([min(timeM_2015(idxs_2015))*.95 240*1.05])
legend('Strategy Stürzebecher 2005','Strategy Cebulla 2006','Strategy Cebulla 2015')
