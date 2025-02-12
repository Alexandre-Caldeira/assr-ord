function CV  = critical_value(alpha,detector,parameters,methods,yin,xin) 


%alpha :n�vel de signific�ncia 
%detector (string): nome do detector 
%parameters: parametros do detector 
% parameters.M -> numero de janelas 
% parameters.L -> banda Lateral
% parameters.N -> n�mero de canais
%methods -> m�todo utilizado para estimar o valor cr�tico
%   pararametrs.
%yin, xin -> sinais utilizadas para as corre��es (cholesky corre��o)
            
%exemplo 
% methods = 'time_Cholesky_corrected'; 
% parameters.M  = 10; 
% y = randn(100,2)*[1 .2; .3, 1];
% critical_value(0.05,'aMSC',parameters,methods,yin) 
%parameters.mistura =='PCA'

%Conferir se par�metros est�o definidos ----------------------
list_detector = {'MSC', 'CSM', 'LFT', 'MMSC', 'aMSC', 'pMSC', 'MCSM', 'aCSM', 'pCSM', 'MLFT', 'aLFT',...
'pLFT', 'MGBT', 'aGBT', 'pGBT','daMSC','dpMSC','daCSM','dpCSM'};
list_methods = {'theoretical', 'Monte_Carlo_default','time_Cholesky_corrected','frequency_Cholesky_corrected'};


if (sum(strcmp(detector,list_detector)) == 0),
    error('Detector not defined'); end

if (sum(strcmp(methods,list_methods)) == 0),
    error('Method not defined'); end

if (alpha<0) || (alpha>1),
    error('The alpha value between 0 and 1'); end
%---------------------------------------------------------------------



%Casos especiais de N=1

if strcmp(detector,'MSC') detector ='aMSC';  parameters.N=1; end
if strcmp(detector,'CSM') detector ='aCSM';  parameters.N=1; end
if strcmp(detector,'LFT') detector ='aLFT';  parameters.N=1; end


if  parameters.N==1
    if strcmp(detector,'daCSM') detector ='aCSM'; end
    if strcmp(detector,'daMSC') detector ='aMSC'; end
    if strcmp(detector,'daLFT') detector ='aLFT'; end
end


%% 1 - Methods:  theoretical 
%Valores te�ricos 
if strcmp(methods,'theoretical') == 1
     
    
    if (strcmp(detector(2:4),'MSC') & (parameters.N==1)) 
    CV = 1 - alpha.^(1./(parameters.M-1)); end

   if (strcmp(detector,'MMSC')) 
    CV = betainv(1-alpha,parameters.N,parameters.M-parameters.N); end
    
  if (strcmp(detector,'MCSM')) 
    CV = chi2inv(1-alpha,2)/(2*parameters.M); end

  if (strcmp(detector,'aCSM')) 
    CV = chi2inv(1-alpha,2)/(2*parameters.M*parameters.N); end
  
  if (strcmp(detector,'MLFT')) 
    CV = finv((1-alpha),2*parameters.N,2*parameters.N*parameters.L); end

end



%% 2- Monte Carlo Default
if strcmp(methods,'Monte_Carlo_default') == 1
    
    %parametros pr�-defindo 
    parameters.fo =8;
    parameters.tj = 32; 
    N = parameters.N;
    
    if(isfield(parameters, 'Nruns'))
        Nruns = parameters.Nruns;
    else
        Nruns = 10000;
    end   
      
    if(strcmp(detector(2),'L')) 
        parameters.M = 1; parameters.tj = 256; parameters.fo = 8*4; parameters.fs = parameters.tj; end
    
    
    % Monte Carlo
    length_signal =  parameters.tj*parameters.M; 
    DV = zeros(Nruns,1); 
    
    for ii = 1:Nruns
       
        if strcmp(detector(3),'B') 
            y = randn(length_signal,N); 
            x = randn(length_signal,N); 
            [value] = mord(detector,y,parameters,x);
        else
            y = randn(length_signal,N); 
            [value] = mord(detector,y,parameters);
        end
        
        
        if strcmp(detector(2),'L') 
            DV(ii)  = value;
        else
            DV(ii) = value(parameters.fo); 
        end 
        
    end
    
    CV = quantile(DV,1-alpha);
      
end



%% Time_Cholesky_corrected
if strcmp(methods,'time_Cholesky_corrected') == 1
    
    %parametrs
    parameters.fo =8;
    parameters.tj = 32; 
    N = parameters.N;    
    if (N ~= size(yin,2)), warning('N ~= number of signals'); end;
    
    if(isfield(parameters, 'Nruns'))
        Nruns = parameters.Nruns;
    else
        Nruns = 10000;
    end   
    
    if(strcmp(detector(2),'L')) 
       parameters.M = 1; parameters.tj = 256; parameters.fo = 8*4; parameters.fs = parameters.tj; end
 
    
    %time_Cholesky_corrected
    if(strcmp(parameters.mistura,'fixa'))
 
        A1 = yin;   
    else
        ymean = mean(yin);
        [nl,nc] = size(yin);
        for j = 1:nc
            yin(:,j) = yin(:,j) - ymean(j);
        end

        if (isfield(parameters, 'mistura')==0 || strcmp(parameters.mistura,'cholesky'))
            Sigma = corr(yin);
            A1 = chol(Sigma); %mixture
        else 
            %parameters.mistura =='PCA'
            [A1] = pca(yin);
            A1 = A1';    
        end
    end
    
    if strcmp(detector(3),'B') 
        xmean = mean(xin);
        [nl,nc] = size(xin);
        for j = 1:nc
        xin(:,j) = xin(:,j) - xmean(j);
        end

        Sigma = corr(xin);
        A2 = chol(Sigma); %mixture 
    end
    
    % Monte Carlo com a matriz de mistura
    length_signal =  parameters.tj*parameters.M; 
    DV = zeros(Nruns,1); 
    
    for ii = 1:Nruns
       
        if strcmp(detector(3),'B') 
            y = randn(length_signal,N)*A1; 
            x = randn(length_signal,N)*A2; 
            [value] = mord(detector,y,parameters,x);
        else
            y = randn(length_signal,N)*A1; 
            y = detrend(y,'constant'); 
            [value] = mord(detector,y,parameters);
        end
        
        
        if strcmp(detector(2),'L') 
            DV(ii)  = value;
        else
            DV(ii) = value(parameters.fo); 
        end 
        
    end %end for
    
    CV = quantile(DV,1-alpha);
   
end



%% Frequency_Cholesky_corrected
if strcmp(methods,'frequency_Cholesky_corrected') == 1
    
    
    %parametros da simula��o
    if (parameters.N ~= size(yin,2)), warning('N ~= number of signals'); end;
    
    if strcmp(detector(2),'L')
        vector_LFT = [-parameters.L/2:parameters.L/2]; 
    end
    %vetor_teste_R = [-9:-8,-4:-1,0,1:4,8:9];
    
    
    if(isfield(parameters, 'Nruns'))
        Nruns = parameters.Nruns;
    else
        Nruns = 10000;
    end   
    
    % identificar o bins correspondente a frequ�ncia analisada 
    parameters.bin =  parameters.tj./parameters.fs*parameters.fo+1;
    
    % estimar o n�mero de janelas
    parameters.M  = floor(size(yin,1)/parameters.tj);
    
    % matriz espectral-cruzada 
    A1freq = Chol_f(yin,parameters.tj);
    
    %Teste-F ????
%     if(strcmp(detector(2),'L')) 
%     parameters.M = 1; parameters.tj = 256; parameters.fo = 8*4; parameters.fs = parameters.tj; end
%  
%     if strcmp(detector(2),'L') %normalized matrix        
%         A1freq = Chol_f_Norm(yin,parameters.tj);
%     end


    if strcmp(detector(3),'B') %normalized matrix        
        A1freq = Chol_f_Norm(yin,parameters.tj);
        A2freq = Chol_f_Norm(xin,parameters.tj);
    end
    
    %Monte Carlo simulation 
    DV = zeros(Nruns,1); 
    for ii = 1:Nruns
        
        if strcmp(detector(3),'B') 
            DV(ii,1) = ORD_freq(A1freq(:,:,parameters.bin),parameters.M,...
                detector,A2freq(:,:,parameters.bin));
        elseif strcmp(detector(2),'L') 
            DV(ii,1) = ORD_freq(A1freq(:,:,parameters.bin+vector_LFT),...
                parameters.M,detector);
        else
            DV(ii,1) = ORD_freq(A1freq(:,:,parameters.bin),parameters.M,...
               detector);  
        end
    end
    CV = quantile(DV,1-alpha);
    
end
