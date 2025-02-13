%clear all,close all,clc
%Analisando ndc sem correção e ndc corrigido
%% analisando a pd
%essa variavel pd corrigido corresponde aos valores da taxa de detecção
%utilizando a correção, foi salvo por colunas as intensidades de 50, 40 e
%30db SPL
%pd_corrigido=TXD; %o txd é obtido do código 'sinal eeg ndc' aplicando ao banco de dados 
%pd_nao_corrigido=TXD; %novamante foi salvo a taxa de detecção sem utilizar a correção em tres colunas
%50,40 e 30db SPL

aumentapd=zeros(size(pd_corrigido,1),size(pd_corrigido,2));
%lembrar de trocar ff para 3
for ff=1:3 %estimar para as tres intensidades
    aumento(:,ff)=pd_corrigido(:,ff)-pd_nao_corrigido(:,ff);
    for ii=1:size(pd_corrigido,1)
        if pd_nao_corrigido(ii,ff)~=0  %para 30db pode ter txd igual a zero 
            aumentapd(ii,ff)=aumento(ii,ff)/pd_nao_corrigido(ii,ff);
       end
    end
end
aumentapd=aumentapd*100;
% por coluna 50db 40db e 30db
%menores intensidades a correção obtem melhor resultado
mean(aumentapd)
%a medinana garante que pelo menos 50% dos resultados estão compreendidos nesses valores 
median(aumentapd)

%em nenhum caso houve piora em relação a utilização da correção
min(aumentapd) 

% novamente 50 40 e 30db 
%distribuição do aumento nas taxas de detecção
figure
hist(aumentapd(:,1),100)
title('Distribuição do aumento da pd 50db')

figure
hist(aumentapd(:,2),100)
title('Distribuição do aumento da pd 40db')

figure
hist(aumentapd(:,3),100)
title('Distribuição do aumento da pd 30db')


figure
boxplot(aumentapd(:,1))
%% analisando o tempo
% A analise do tempo mostrou que não teve grandes diferenças utilizando a correção
%sempre bom frisar que em nenhum caso piorou em relação a não correção

%tempo_NDC_semCorrecao=timeMedio; %novamente a varivel timeMedio pode ser obtido aplicando a estrategia
%ao banco de dados
%tempo_NDC_alfaCorrigdo=timeMedio;

reduztempo=zeros(size(tempo_NDC_alfaCorrigdo,1),size(tempo_NDC_semCorrecao,2));
for gg=1:3  %para as tres insidades
    reducao(:,gg)=tempo_NDC_alfaCorrigdo(:,gg)-tempo_NDC_semCorrecao(:,gg); 
   for hh=1:size(tempo_NDC_semCorrecao,1)
        reduztempo(hh,gg)=reducao(hh,gg)/tempo_NDC_semCorrecao(hh,gg);
   end
end
%valores negativos indicam a correção mais rápida
%por ordem 50db 40db e 30db
reduztempo=reduztempo*100;
mean(reduztempo)
median(reduztempo)

max(reduztempo)  %em nenhum caso teve aumento no tempo de exame

%distribuição da redução do tempo
figure
hist(reduztempo(:,1),100)
title('Distribuição da redução do tempo 50db')

figure
hist(reduztempo(:,2),100)
title('Distribuição da redução do tempo 40db')

figure
hist(reduztempo(:,3),100)
title('Distribuição da redução do tempo 30db')

%% plotando as fronteiras pareto
%criando a figura

TXD=auxpd50db(:,1); 
timeM=auxtime50db(:,1);
%essa varivel 'auxpd50db' é uma varivel que contem as taxas de detecção
%para todos os conjuntos de parametros, ela é dividida em duas colunas, na
%primeira coluna são as taxas de detecção utilizando a correção e na
%segunda coluna é a taxa de detecção sem a correção a mesma análise vale
%para auxtime50db, nese caso precisa de rodar novremnte utilizando as
%intendsidades de 40 e 30db

figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');  %estudar
plot([0 1]*100, [Mmax Mmax],'-.k','linewidth',1) %linha do mmax 240
plot([TXD(1) TXD(1)], [min(timeM) max(timeM)],'-.k','linewidth',1.1) %linha do teste unico
%Teve alteração nessa parte

for hh=1:2 %primeira rodada plota ndc com correção, na segunda sem correção
    TXD=auxpd50db(:,hh);
    timeM=auxtime50db(:,hh);
    
    %plota os segmentos em pontos
     if hh==1   
        for ii = 1:size(parametros,1)
            plot(TXD(ii),timeM(ii),'.k','Markersize',4,'DisplayName',[num2str(parametros(ii,1)) '-' num2str(parametros(ii,2))])
        end
     end
      %if hh==2 %caso queira plotar os protocolos do ndc sem correção   
       % for ii = 1:size(parametros,1)
       %     plot(TXD(ii),timeM(ii),'Marker','.','Color',[0 0.4470 0.7410]','Markersize',6,'DisplayName',[num2str(parametros(ii,1)) '-' num2str(parametros(ii,2))])
       % end
      %end
    [ p, idxs] = paretoFront([TXD,(-timeM)] );
    auxL = p(:,1)<0.5; 
    p(auxL,:) = [];
    idxs(auxL,:) = [];

    [~,ind] = sort(p(:,1));
    p = p(ind,:);
    idxs = idxs(ind,:);

    if hh==1
        plot(TXD(idxs),timeM(idxs),'-ok','Markersize',8,'linewidth',1.8) %plot da fronteira
        %legend('NDC Com Alfa Corrigido');
    else
        hold on
        plot(TXD(idxs),timeM(idxs),'-^r','Markersize',8,'linewidth',1.8) %plot da fronteira
        %legend('NDC Sem Correção');
    end
    for ii = 1:size(idxs,1)
     [I] = find((TXD == TXD(idxs(ii))) & (timeM==timeM(idxs(ii))));
     fprintf('\nPD = %f Tempo = %f ',TXD(idxs(ii)),timeM(idxs(ii))); 
     fprintf(' NI = %d ',length(I)); 
     I = I(1); 
        for jj = I  
        
            fprintf(' - Buffer:%d, M_step:%d', parametros(jj,1),parametros(jj,2)); 
            if hh==1
                text(TXD(jj),timeM(idxs(ii))*.995,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) '\}' ],'Fontsize',9); %plotar o texto dos procolos
            else
                 text(TXD(jj)*.979,timeM(idxs(ii))*1.01,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) '\}' ],'Fontsize',9); %plotar o texto dos procolos
            end
        end    
    end
end
hold off
xlabel('Taxa de Detecção(%)','fontsize',12); 
ylabel('Tempo de Exame Médio (Segundos)','fontsize',12);
text(max(auxpd50db(:,1)),245,'NDC com ajuste do nível de significância','Color','black','FontSize',10)
text(max(auxpd50db(:,1)),248,'NDC sem correção','Color','red','FontSize',10)
fprintf('\n');  
xlim([min(TXD(idxs))*.95,max(TXD(idxs))*1.05])
ylim([min(timeM(idxs))*.95 Mmax*1.05])
title('Intensidade de 50 db SPL','Fontsize',12)

%% Plotando os protocolos que obtem o melhor desempenho

x=parametros(:,1:2); %conjunto Mmin,Mstep,Mmax
x(:,3)=aumentapd(:,1);%50db
x(:,4)=aumentapd(:,2);%40db
x(:,5)=aumentapd(:,3);%30db

for hh=3:5
    for ii = 1:size(x,1)
    
        if x(ii,hh)<10  %menor que 10% de aumento
            cor(ii,:) = [0 0.4470 0.7410];
            tamanho(ii,1) = 10;
        elseif x(ii,hh)<25
            cor(ii,:) = [0.4660 0.6740 0.1880];
            tamanho(ii,1) = 40;
        elseif  x(ii,hh)<50
            cor(ii,:) = [0.9290 0.6940 0.1250];
            tamanho(ii,1) = 55;
        else
            cor(ii,:) = [0.6350 0.0780 0.1840];
            tamanho(ii,1) = 75;
        end
        
    end
figure,
s=tamanho;
scatter3(x(:,1),x(:,2),x(:,hh),s,cor)
xlabel('M_{MIN}')
ylabel('M_{STEP}')
zlabel('Aumento PD(%)')
%legend({'Line 1','Line 2','Line 3'})
view(0,90)    
    
   
end
%a vista superior fornece a imagem semelhenate ao exemplo proposto
figure,
s=tamanho;
scatter3(x(:,1),x(:,2),x(:,3),s,cor)
xlabel('M_{MIN}')
ylabel('M_{STEP}')
zlabel('Aumento PD(%)')
%legend({'Line 1','Line 2','Line 3'})
title('Protocolos de detecção 50db')

%colorbar

%% salvando variaveis
%save(['NDC_AlfaCorrigido_CSM_Mmax' num2str(Mmax) '_alfa_'  num2str(alfa) '_FPdesejado' num2str(FP_desejado) '.mat'],'alfa_corrigido', ...
%    'NDC_minimo','P', 'nRuns','alfa','FP_desejado')


%informacao={'NDC_alfa_corrigido', 'NDC_sem_correção','Variaveis_para_plotar_fronteira_pareto'};
%save(['timeM_30db_ndc_sem_correcao''.mat'],'timeM')
%save(['pareto_font_NdcAlfaCorrigido_NdcSemCorrecao''.mat'],'auxpd30db','auxpd40db','auxpd50db','auxtime30db','auxtime40db','auxtime50db','Mmax','parametros');
