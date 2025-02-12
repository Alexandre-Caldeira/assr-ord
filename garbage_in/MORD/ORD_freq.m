function ORD = ORD_freq(H,M,s,Hr)
%gerar fun��o valores cr�ticos 

%parametros 
%Nome da fun��o
%M n�mero da janela
%N ser� obtido pela matriz de mistura na frequ�ncia. 
%Nruns = 1; %sempre igual a 1. 
% s = ['aMSC'];
% M = 30; 
% H = 1; 

N = size(H,1);
%gear valores do detector para a matriz na frequ�ncia
if strcmp(s,'aMSC') == 1 
    
    yfft = (randn(M,N)+j*randn(M,N))*H;
    MSC1 = zeros(1,N);
    for jj = 1:N

        MSC1(1,jj) = abs(sum(yfft(1:M,jj),1)).^2./(M*sum(abs(yfft(1:M,jj)).^2,1));

    end
    ORD = mean(MSC1,2); %aMSC
end

if strcmp(s,'pMSC') == 1 
    
    yfft = (randn(M,N)+j*randn(M,N))*H;
    MSC1 = zeros(1,N);
    for jj = 1:N

        MSC1(1,jj) = abs(sum(yfft(1:M,jj),1)).^2./(M*sum(abs(yfft(1:M,jj)).^2,1));

    end
    ORD = (prod(MSC1,2)).^(1/N); %aMSC
end


if strcmp(s,'pCSM') == 1 
    
    yfft = (randn(M,N)+j*randn(M,N))*H;
    yfft = yfft./abs(yfft);
    CSM1 = zeros(1,N);
    for jj = 1:N

        CSM1(1,jj) = abs(sum(yfft(1:M,jj),1)).^2./(M*sum(abs(yfft(1:M,jj)).^2,1));

    end
    ORD = (prod(CSM1,2)).^(1/N); %aMSC
end


if strcmp(s,'aCSM') == 1 
    
    yfft = (randn(M,N)+j*randn(M,N))*H;
    yfft = yfft./abs(yfft);
    CSM1 = zeros(1,N);
    for jj = 1:N

        CSM1(1,jj) = abs(sum(yfft(1:M,jj),1)).^2./(M*sum(abs(yfft(1:M,jj)).^2,1));

    end
    ORD = mean(CSM1,2); %aMSC
end


if strcmp(s,'MCSM') == 1 
    yfft = (randn(M,N)+j*randn(M,N))*H;
    teta = angle(yfft);
    C = cos(teta);
    S = sin(teta);
    Cmed = mean(C,2);
    Smed = mean(S,2);
    temp1 = atan( (Smed.*(Cmed<0))./Cmed)+pi*(Cmed<0); %calculates mean teta matrix only for Cmed<0, zero-padding other positions
    temp2 = atan( (Smed.*(Cmed>=0))./Cmed); %calculates mean teta matrix only for Cmed>0, zero-padding other positions
    teta_med = temp1 + temp2; 
    ORD = (1/M^2)*sum(cos(teta_med)).^2+(1/M^2)*sum(sin(teta_med)).^2; %multiple CSM
end



if strcmp(s,'aLFT') == 1 
    
    L = size(H,3);
    yfft = zeros(L,N);
    for ii = 1:L 
        yfft(ii,1:N) = (randn(1,N)+j*randn(1,N))*H(:,:,ii);
    end
    
    pfo = (L-1)/2+1;
    L = L-1;
    Y = abs(yfft); 
    Yfo = Y(pfo,:); 
    Yfn = Y;
    Yfn(pfo,:) = []; 

    % Compute F value
    ORD = mean((Yfo.^2)./(1/L*sum(Yfn.^2,1)));
    
end


if strcmp(s,'pLFT') == 1 
    
    L = size(H,3);
    yfft = zeros(L,N);
    for ii = 1:L 
        yfft(ii,1:N) = (randn(1,N)+j*randn(1,N))*H(:,:,ii);
    end
    
    pfo = (L-1)/2+1;
    L = L-1;
    Y = abs(yfft); 
    Yfo = Y(pfo,:); 
    Yfn = Y;
    Yfn(pfo,:) = []; 

    % Compute F value
    ORD = (prod((Yfo.^2)./(1/L*sum(Yfn.^2,1)))).^(1/N);
    
end

if strcmp(s,'MLFT') == 1 
    
    L = size(H,3);
    yfft = zeros(L,N);
    for ii = 1:L 
        yfft(ii,1:N) = (randn(1,N)+j*randn(1,N))*H(:,:,ii);
    end
    
    pfo = (L-1)/2+1;
    L = L-1;
    Y = abs(yfft); 
    Yfo = Y(pfo,:); 
    Yfn = Y;
    Yfn(pfo,:) = []; 

    
    ORD = sum(Yfo.^2)/(sum(1/L*sum(Yfn.^2,1)));
    
end


if strcmp(s,'aGBT') == 1 
    
    
    Svv = sum((abs((randn(M,N)+j*randn(M,N))*H)).^2,1);
    Snm = sum((abs((randn(M,N)+j*randn(M,N))*Hr)).^2,1);

    ORD  = mean(Svv./(Svv+Snm),2);

end

if strcmp(s,'pGBT') == 1 
    
    
    Svv = sum((abs((randn(M,N)+j*randn(M,N))*H)).^2,1);
    Snm = sum((abs((randn(M,N)+j*randn(M,N))*Hr)).^2,1);

    ORD  = (prod(Svv./(Svv+Snm),2)).^(1/N);

end


if strcmp(s,'MGBT') == 1 
    
    
    Svv = sum(sum((abs((randn(M,N)+j*randn(M,N))*H)).^2,1));
    Snm = sum(sum((abs((randn(M,N)+j*randn(M,N))*Hr)).^2,1));

    ORD  = Svv./(Svv+Snm);

end


