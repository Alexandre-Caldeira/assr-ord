%Código que Avalia o ruído do canais 


clear all
close all
clc 


%vetor dos voluntários 
Vvoluntario = {'Abdon';'Ana';'BBB';'Colatina';'Erick';'Luciana';...
    'Sombra';'Quenaz';'Vinicius';'Sacola';'Wreikson'}; %vetor dos voluntário 

%vetor da intensidade 
Vintensidade = {'70';'60';'50';'40';'30'}; 
load('eletrodos.mat')

%% --------------------------------------------------
remoc = [.1]/200; 
sweep = 1; %número de janelas em um sweep 
freqM = 150; 
fcInferior = 78;
fcSuperior = 98;
lw = 1;

%% 
A = nan(8,1000,5,16,11); %vetor amplitude 
R = nan(freqM,1000,16,11);
% salvar as variáveis 
% 1 trechos
% 2  Moduladora 
% 3 Combinação de eletrodos  
cont = 1; 
for cont_vol = 1:size(Vvoluntario,1)

    voluntario = cell2mat(Vvoluntario(cont_vol,:));
    
    for cont_int = 1:size(Vintensidade,1) 
        
        intensidade = cell2mat(Vintensidade(cont_int,:));
        load([voluntario '_'  intensidade 'dB'], 'x','Fs','binsM','freqEstim')   
  
        nfft = Fs;%1segundo de sinal 
         
        %retirar componente DC por janela (fiz isso pq no processamento em
        %tempo real é por janela)
        x = x - repmat(mean(x),nfft,1); %tirar a média de cada de cada trecho - devido a remoção
        %encontrar o valor máximo por canal 
        Vmax = squeeze(max(max(abs(x)),[],3));
        ind =[];
        ind = Vmax>remoc;
        x(:,ind,:) = []; 
        
        
        [b,a] = butter(4,[fcInferior/(Fs/2), fcSuperior/(Fs/2)]);
        x = reshape(x(:,1:(floor(size(x,2)/sweep)*sweep),:),[],16);
        x = filter(b,a,x); 
        x = reshape(x,nfft*sweep,[],16);
        
        %scluir os dois primeiros segundos do inicio da coleta 
        x(:,1,:) =[]; 
        
        xmedia = squeeze(mean(x,2));
        xx = fft(x)./(size(x(:,cont),1)); 
         A(:,1:size(xx,2),cont_int,:,cont_vol) = xx((binsM-1)*sweep+1,:,:); 
                 
 
    end
 %   pause 
 %   close all
end


%% Fazer com EEG espontaneo 

Vvoluntario = {'AbdonESP';'AnaESP';'BBBESP';'ColatinaESP';'ErickESP';'LucianaESP';...
    'SombraESP';'QuenazESP';'ViniciusESP';'SacolaESP';'WreiksonESP'}; %vetor dos voluntáriovetor dos voluntário 

cont = 1;   
for cont_vol = 1:size(Vvoluntario,1)

       voluntario = cell2mat(Vvoluntario(cont_vol,:));

       load([voluntario], 'x','Fs','binsM','freqEstim')   

             nfft = Fs;%1segundo de sinal 
         
        %retirar componente DC por janela (fiz isso pq no processamento em
        %tempo real é por janela)
        x = x - repmat(mean(x),nfft,1); %tirar a média de cada de cada trecho - devido a remoção
        %encontrar o valor máximo por canal 
        Vmax = squeeze(max(max(abs(x)),[],3));
        ind =[];
        ind = Vmax>remoc;
        x(:,ind,:) = []; 
        
        
        [b,a] = butter(4,[fcInferior/(Fs/2), fcSuperior/(Fs/2)]);
        x = reshape(x(:,1:(floor(size(x,2)/sweep)*sweep),:),[],16);
        x = filter(b,a,x); 
        x = reshape(x,nfft*sweep,[],16);
        
        %scluir os dois primeiros segundos do inicio da coleta 
        x(:,1,:) =[]; 
        
        xmedia = squeeze(mean(x,2));
        xx = fft(x)./(size(x(:,cont),1)); 
        R(:,1:size(xx,2),:,cont_vol) = xx(1:freqM,:,:); 
    
        
end    

R = R.*1e9;
A = A.*1e9;
%% plotar a figuras do ruído 
%Msweep deve ser igual a 1. 
close all
vol = [1,2,3,4,5,8,10,11];
% vol = [2];

%figure 1 
%plotar o ruído em função da frequêncai 
Mint = [5,10,40,70,100]; 
mm = ['kp';'kv';'ks';'ko';'k<']
figure('Name', 'Ruído em função da moduladora');
freq= 0:freqM-1
lw =1.5;
mk = 8;


ll = [ repmat('M=',size(Mint,2),1),num2str(Mint')];

for ii = 1:size(Mint,2)
    
   hold on 
    plot(binsM-1,mean(mean(abs(mean(R(binsM,1:Mint(ii),:,vol),2)),3),4),mm(ii,:),...
        'markersize',mk,'MarkerFaceColor','k');
   
end

legend({ll},'Location','eastoutside','EdgeColor','none','fontsize',12)
for ii = 1:size(Mint,2)
   hold on 
   plot(freq,mean(mean(abs(mean(R(:,1:Mint(ii),:,vol),2)),3),4),'k','LineWidth',lw)   
end


xlabel('Frequência (Hz)','fontsize',12)
ylabel('Amplitude (nV)','fontsize',12)
set(gca,'YMinorTick','on')
box off
xlim([79 98])


% R = nan(freqM,1000,16,11);



%% 
vol = [1,2,3,4,5,8,11,10]; 
vol = [1,2,3,4,5,8,10,11]
% vol = [5]

figure('Name','Ruído em função janela')
lw = 1.5;
canal = 1;
Mint = 100;
rmedia=[];
for ii = 1:Mint(end)
   
    rmedia(ii) = mean(mean(mean(abs(mean(R(binsM,1:ii,canal,vol),2)),3),4),1);
end

plot(1:Mint,rmedia,'k','LineWidth',lw)
xlim([1 Mint])
xlabel('M','fontsize',12)
ylabel('Amplitude (nV)','fontsize',12)
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on')
box off

A50=[];
for ii = 1:Mint(end)
   
    A50(ii) = mean(mean(mean(abs(mean(A(:,1:ii,3,canal ,vol),2)),4),5),1);
end

A40=[];
for ii = 1:Mint(end)
   
    A40(ii) = mean(mean(mean(abs(mean(A(:,1:ii,4,canal ,vol),2)),4),5),1);
end

A30=[];
for ii = 1:Mint(end)
   
    A30(ii) = mean(mean(mean(abs(mean(A(:,1:ii,5,canal ,vol),2)),4),5),1);
end



figure 
plot(1:Mint(end),rmedia,'k','LineWidth',2*lw)
xlim([1 Mint])
hold on 
plot(1:Mint(end),A50,'k:','LineWidth',lw)
plot(1:Mint(end),A40,'k-.','LineWidth',lw)
%plot(1:Mint(end),A30,'r','LineWidth',lw)

xlabel('M','fontsize',12)
ylabel('Amplitude (nV)','fontsize',12)
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on')
box off
legend({'Espontâneo';'50dB';'40dB'}, 'Location','best',...
    'edgecolor','none','fontsize',12)


%% 
MM = [50,300,250,420,420]; 
MM = [50,300,250,100,100]; 
MM = [100,120,100,100,100];
canal = 1; 
%vol = [5]
% vol =1:11; 
vol = [1,2,3,4,5,8,10,11];
% vol =2;
lw = 1.5
for ii = 3:size(MM,2)
    AA(:,ii) = mean(mean(mean(abs(mean(A(:,1:MM(ii),ii,canal ,vol),2)),4),5),3);
end


figure 
plot(freqEstim,AA(:,3:end),'LineWidth',lw)
xlabel('Frequência (Hz)','fontsize',12)
ylabel('Amplitude (nV)','fontsize',12)
set(gca,'YMinorTick','on')
box off
% legend({'70','60','50','40','30'})
legend({'70','60','50','40','30'})

% xlim([79 98])
hold on
plot(binsM-1,mean(mean(abs(mean(R(binsM,1:100,canal,vol),2)),3),4),'-k','LineWidth',lw);
legend({'50 dB','40 dB','30 dB','Esp'},'Location','eastoutside','EdgeColor','none','fontsize',12)
title(['M=100'],'fontsize',12)
 

%% 
AA =[];
MM = [100,100,100,100,100];
canal = 14; 
%vol = [5]
% vol =1:11; 
vol = [1,2,3,4,5,8,10,11];
% vol =2;
lw = 2
for ii = 3:size(MM,2)
    AA(:,ii) = mean(mean(mean(abs(mean(A(:,1:MM(ii),ii,canal ,vol),2)),4),5),3);
end


figure 
plot(freqEstim,AA(:,3:end),'LineWidth',lw)
xlabel('Frequência (Hz)','fontsize',12)
ylabel('Amplitude (nV)','fontsize',12)
set(gca,'YMinorTick','on')
box off
% legend({'70','60','50','40','30'})
% legend({ll},'Location','eastoutside','EdgeColor','none','fontsize',12)



MM = [100,100,250,400,400];
lw = 1.5
for ii = 3:size(MM,2)
    AA(:,ii) = mean(mean(mean(abs(mean(A(:,1:MM(ii),ii,canal ,vol),2)),4),5),3);
end
hold on 
plot(freqEstim,AA(:,3:end),':','LineWidth',lw)


% xlim([79 98])
% hold on
% plot(binsM-1,mean(mean(abs(mean(R(binsM,1:100,canal,vol),2)),3),4),'-k');
%  



%% 
canal = 16;
Mint = 100;
figure('name',['Canal=' num2str(canal) ' M=' num2str(Mint)])


ll = {'500 Hz';'1000 Hz';'2000 Hz';'4000 Hz'}; 
for kk =1:4 
    %vol = [1,2,3,4,5,8,11,10]; 
    vol = [1,2,3,4,5,8,10,11];
    % vol = [5]

%     figure('Name','Ruído em função janela')
    lw = 1.5;
    rmedia=[];
    for ii = 1:Mint(end)
        rmedia(ii) = mean(mean(mean(abs(mean(R(binsM((2*kk-1):(2*kk)),1:ii,canal,vol),2)),3),4),1);
    end

%     plot(1:Mint,rmedia,'k','LineWidth',lw)
%     xlim([1 Mint])
%     xlabel('M','fontsize',12)
%     ylabel('Amplitude (nV)','fontsize',12)
%     set(gca,'YMinorTick','on')
%     set(gca,'XMinorTick','on')
%     box off

    A50=[];
    for ii = 1:Mint(end)

        A50(ii) = mean(mean(mean(abs(mean(A((2*kk-1):(2*kk),1:ii,3,canal ,vol),2)),4),5),1);
    end

    A40=[];
    for ii = 1:Mint(end)

        A40(ii) = mean(mean(mean(abs(mean(A((2*kk-1):(2*kk),1:ii,4,canal ,vol),2)),4),5),1);
    end

    A30=[];
    for ii = 1:Mint(end)

        A30(ii) = mean(mean(mean(abs(mean(A((2*kk-1):(2*kk),1:ii,5,canal ,vol),2)),4),5),1);
    end



    subplot(2,2,kk)
 
    plot(1:Mint(end),rmedia,'k','LineWidth',2*lw)
    xlim([1 Mint])
       ylim([0 250])
    hold on 
    plot(1:Mint(end),A50,'k:','LineWidth',lw)
    plot(1:Mint(end),A40,'k-.','LineWidth',lw)
%     plot(1:Mint(end),A30,'r','LineWidth',lw)

    title(ll(kk),'fontsize',12)
    xlabel('M','fontsize',12)
    ylabel('Amplitude (nV)','fontsize',12)
    set(gca,'YMinorTick','on')
    set(gca,'XMinorTick','on')
    box off
    legend({'Espontâneo';'50dB';'40dB'}, 'Location','NorthEast' ,...
        'edgecolor','none','fontsize',10)
end


%% 
canal = 16;
Mint = 100;
figure('name',['Canal=' num2str(canal) ' M=' num2str(Mint)])


ll = {'500 Hz';'1000 Hz';'2000 Hz';'4000 Hz'}; 
for kk =1:4 
    %vol = [1,2,3,4,5,8,11,10]; 
    vol = [1,2,3,4,5,8,10,11];
    % vol = [5]

%     figure('Name','Ruído em função janela')
    lw = 1.5;
    rmedia=[];
    for ii = 1:Mint(end)
        rmedia(ii) = mean(mean(mean(abs(mean(R(binsM((2*kk-1):(2*kk)),1:ii,canal,vol),2)),3),4),1);
    end

%     plot(1:Mint,rmedia,'k','LineWidth',lw)
%     xlim([1 Mint])
%     xlabel('M','fontsize',12)
%     ylabel('Amplitude (nV)','fontsize',12)
%     set(gca,'YMinorTick','on')
%     set(gca,'XMinorTick','on')
%     box off

    A50=[];
    for ii = 1:Mint(end)

        A50(ii) = mean(mean(mean(abs(mean(A((2*kk-1):(2*kk),1:ii,3,canal ,vol),2)),4),5),1);
    end

    A40=[];
    for ii = 1:Mint(end)

        A40(ii) = mean(mean(mean(abs(mean(A((2*kk-1):(2*kk),1:ii,4,canal ,vol),2)),4),5),1);
    end

    A30=[];
    for ii = 1:Mint(end)

        A30(ii) = mean(mean(mean(abs(mean(A((2*kk-1):(2*kk),1:ii,5,canal ,vol),2)),4),5),1);
    end



    subplot(2,2,kk)
 
    plot(1:Mint(end),rmedia,'k','LineWidth',2*lw)
    xlim([1 Mint])
       ylim([0 250])
    hold on 
    plot(1:Mint(end),A50,'k:','LineWidth',lw)
%    plot(1:Mint(end),A40,'k-.','LineWidth',lw)
%     plot(1:Mint(end),A30,'r','LineWidth',lw)

    title(ll(kk),'fontsize',12)
    xlabel('M','fontsize',12)
    ylabel('Amplitude (nV)','fontsize',12)
    set(gca,'YMinorTick','on')
    set(gca,'XMinorTick','on')
    box off
    legend({'Espontâneo';'50dB'}, 'Location','NorthEast' ,...
        'edgecolor','none','fontsize',10)
end


