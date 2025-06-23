function PDtime(detector,methods,Mmax,A) 





%% Parametros do protocolo de detecção. 

alpha = 0.05;
parameters.Nruns = 100000; %valor crítico 
parameters.M = Mmax;
parameters.L = 12;
N = size(A,1);
parameters.N = N;
Nint = 30000; % para estimar a taxa de PD
vetor_SNR = [-60:2:10]; 

%parametros filtro 
a_filter = 1;
b_filter = [1 .1 .01];

%%

NFFT = 32; %lowest base two with good results. 
parameters.fs = NFFT; 
parameters.tj = NFFT; 

freq = round(NFFT*.1);
s = [cos(2*pi*freq*[1:(Mmax*NFFT)]/(NFFT))]';
 %   s = cos(2*pi*freq*linspace(0,M,M*NFFT));
sigma_n = 2/NFFT;


%methods = ['time_Cholesky_corrected'];
parameters.mistura = 'fixa';
if (strcmp(detector,'MMSC') ||strcmp(detector,'MCSM'))
    crit = critical_value(alpha,detector,parameters,'theoretical',A); 
else 
    crit = critical_value(alpha,detector,parameters,methods,A); 
end
PD = zeros(size(vetor_SNR,2),1);
for jj =1:size(vetor_SNR,2)
    SNR = 10^(vetor_SNR(jj)/10);
    

    vPD = zeros(Nint,1);
    for ii = 1:Nint %número de repetições
        
           Noise = randn(size(s,1),N);
           %filtro ------------------------------
%             Noise = randn(size(s,1)*2,N);
%            Noise = filter(b_filter,a_filter,Noise);
%            Noise(1:round(size(s,1)/3),:) =[];
%            Noise= Noise(1:(size(s,1)),:);
           %-------------------------------------
           
           CorNoise = Noise*A;                
           y = repmat(sqrt(4*sigma_n*SNR/NFFT)*s,1,N) + sqrt(sigma_n).* CorNoise;

      
          if(strcmp(detector(2),'L')) 
              
                  parameters.fo  = freq; %valor da freq.
                  aux = mord(detector,y,parameters);
                 vPD(ii) = aux>crit;
          else
               aux = mord(detector,y,parameters);
               vPD(ii) = aux(freq+1)>crit;
          end
    end
    
    PD(jj) = mean(vPD);
end
      
save([ detector '_' methods '_' num2str(Mmax)],'PD','vetor_SNR','A','parameters','a_filter','b_filter');
%save([ detector '_' methods '_' num2str(Mmax) '_filter'],'PD','vetor_SNR','A','parameters','a_filter','b_filter');



