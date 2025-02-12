function H = Chol_f(y,L)

% function [x,F] = Chol_f(y,L);
% Esta  rotina retorna uma matriz cujas colunas sao as trasformadas de
% Fourier das misturas dos sinais nas colunas de y de acordo a matriz de
% espectros cruzados Syy.
%
% Parametros de entrada:    y       => matriz cujas colunas sao os sinais considerados
% 			                                
[tamsinal,N] = size(y);
nfft = fix(L/2);     
M = fix(tamsinal/L); 
y = y(1:M*L,:);      %joga fora amostras que est�o sobrando. Evitar janelas com menos pontos
for i = 1:N,
    Y(:,:,i)=fft(reshape(y(:,i),L,M)); %calcula matriz 3-D cuja i-esima fatia e� a DFT das M janelas de yi[k], i=1,2,...N 
end
Y = Y(1:nfft+1,:,:);%s� retorna metade mais um dos valores, pois o restante � complexo conjugado dos primeiros

%Y = Y./repmat((mean(abs(Y),2)),1,M);


for i=1:N, 
   for j=i:N, 
      Syy(i,j,:) = sum( (conj(Y(:,:,i)).*Y(:,:,j)).' );
   end 
end 
for i = 1:N, 
   for j = 1:i-1, 
      Syy(i,j,:) = conj(Syy(j,i,:)); 
   end 
end 

%xr = randn(M,N,nfft);
%xi = randn(M,N,nfft);
%Syy_r = real(Syy);
%Syy_i = imag(Syy);
%x = zeros(M,N,nfft);

H = zeros(N,N,nfft);
for i =2:nfft,    % N�o incluir o dc
    %Cr = chol(Syy_r(:,:,i));
    %Ci = chol(Syy_i(:,:,i));
   % xr(:,:,i) = xr(:,:,i)*Cr;
   % xi(:,:,i) = xi(:,:,i)*Ci;
 %  x(:,:,i) = x(:,:,i) + (randn(M,N) + j*(rand(M,N)))/sqrt(2);
   H(:,:,i) = chol(Syy(:,:,i));
 %  x(:,:,i) = x(:,:,i)*C;
end
%x = xr + j*xi;
