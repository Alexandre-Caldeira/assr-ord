function [betaN,F,betaNcrit] = MGBT(y,x,tj,fs,alpha)
%
% Multivariate Global Beta test
%
%Extension of the global Beta test to multi channel analysis. This is a MORD 
%technique that makes use of only magnitude information and can be used as a 
%detector of hidden periodicities in noise, since N pairs of real-valued signals are
%available. The pure noise and signal+noise time-series form a pair. 
%The method is described in "not yet."
%
%Sintaxes:
%
%betaN = MGBT(y,x,tj) => MGFT spectrum.
%[betaN,F] = MGBT(y,x,tj,fs) => MGBT spectrum and frequency vector. 
%[betaN,F,betaNcrit] = MGBT(y,x,tj,fs,alpha) =>  returns also the theoretical 
%critical value.
%
%Input Parameters:
%
%y => matrix whose columns are the signals.  
%x => matrix whose columns are the noise.  
%tj => number of points of each epoch in which the signals will be divided.
%fs => sample rate of signals.
%alpha => significance level, e.g. alpha = 0.05.
%
%Example: MGBT using two signals 
%
% y1 = awgn(sin(2*pi*21*(linspace(0,10,1000))),-15,'measured','db')';
% y2 = awgn(sin(2*pi*21*(linspace(0,10,1000))),-16,'measured','db')';
% x = randn(1000,2);
% [betaN,F,betaNcrit] = MGBT([y1 y2],x,100,100,0.05);
% figure;plot(F,betaN,'b',21,1.01*betaN(22),'-rv', [0 50],[betaNcrit betaNcrit],'k--')
% axis([0 F(end) 0 1]);xlabel('Frequency (Hz)');ylabel('MGBT')

% 2)falso alarme 
% nfft = 50000; 
% M = 10;
% Mx = 9; 
% y = randn(nfft*M,2);
% x = randn(nfft*Mx,2); 
% [betaN,F,MGFTcrit] = MGFT(y,x,nfft,nfft,0.05);
% figure;plot(F,betaN,'b',[0 F(end)],[MGFTcrit MGFTcrit],'k--')
% axis([0 F(end) 0 1]);xlabel('Frequency(Hz)');ylabel('MGFT')
% text(10,min(MGFTcrit*1.2,.98), ['FA = ' num2str(mean(double(betaN>MGFTcrit))*100) '%'],'fontsize',12)
% fprintf('\nFA = %f\n', mean(double(betaN>MGFTcrit)))


nf = fix(tj/2)+1;
N = size(y,2); %numero de sinais do
Nx = size(x,2); %numero de sinais do
if N ~=Nx  error('different number of signals (N)'); end 

M = fix(size(y,1)/(tj));   %determina numero de segmentos
Mx = fix(size(x,1)/(tj));   %determina numero de segmentos
y = y(1:tj*M,:);
x = x(1:tj*Mx,:);
%normaliar pela energia de cada sinal 
% y = y./repmat(std(y),tj*M,1);
% x = x./repmat(std(x),tj*Mx,1);
y = reshape(y,tj,M,N);
x = reshape(x,tj,Mx,N);
%normalizar os sinais pela energia total dos sinais 
y = y./repmat(std(y),tj,1,1); %normalizar pela janela 
x = x./repmat(std(x),tj,1,1);


Y = abs(fft(y)).^2;
X = abs(fft(x)).^2;
Svv = sum(sum(Y(1:nf,:,:),2),3);
Snm = sum(sum(X(1:nf,:,:),2),3);


betaN(1,:) = Svv./(Svv+Snm); 


if nargin > 3
    F = (0:(nf-1))*fs/tj; %frequency vector
    if nargin == 5
        betaNcrit = betainv(1-alpha,M*N,Mx*N); %critical value
    end
end   




