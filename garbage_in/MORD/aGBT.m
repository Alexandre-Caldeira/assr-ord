function [abeta,F,abetacrit] = aGBT(y,x,tj,fs,alpha)
%
% Additive local Beta test.
%
% Extension of the Beta test to multi channel analysis as the sum
% GBTs of each channel. This is a MORD technique that
% makes use only magnitude information and can be used as a detector of
% hidden periodicities in noise, since N real-valued signals are available.
% The method is described in "not yet"
%
%Sintaxes:
%
%abeta = aGBT(y,x,tj) => aGBT spectrum.
%[abeta,F] = aGFT(y,x,tj,fs) => aGBT spectrum plus frequency vector.
%[abeta,F,abetacrit] = aGBT(y,tj,x,fs,alpha) =>  returns also the simulated
%critical value.
%
%Input Parameters:
%
%y => matrix whose columns are the signals.
%x => matrix whose columns are the noises.
%tj => number of points of each epoch in which the signals will be divided.
%fs => sample rate of signals.
%alpha => significance level, e.g. alpha = 0.05.
%
%Example: aGFT using two signals
%
% y1 = awgn(sin(2*pi*21*(linspace(0,10,1000))),-15,'measured','db')';
% y2 = awgn(sin(2*pi*21*(linspace(0,10,1000))),-16,'measured','db')';
% x = randn(1000,2);
% [abeta,F,abetacrit] = aGBT([y1 y2],x,100,100,0.05);
% figure;plot(F,abeta,'b',21,1.01*abeta(22),'-rv', [0 50],[abetacrit abetacrit],'k--')
% axis([0 F(end) 0 1]);xlabel('Frequency (Hz)');ylabel('aGBT')

% %2) falso alarme
% nfft = 50000;
% M = 10;
% y = randn(nfft*M,2);
% x = randn(nfft*M,2);
% [abetaN,F,aGFTcrit] = aGFT(y,x,nfft,nfft,0.05);
% figure;plot(F,abetaN,'b',[0 F(end)],[aGFTcrit aGFTcrit],'k--')
% axis([0 F(end) 0 1]);xlabel('Frequency(Hz)');ylabel('MGFT')
% text(10,min(aGFTcrit*1.2,.98), ['FA = ' num2str(mean(double(abetaN>aGFTcrit))*100) '%'],'fontsize',12)
% fprintf('\nFA = %f\n', mean(double(abetaN>aGFTcrit)))


% 3) Falso alarme para diferentes graus de liberdade
% nfft = 100000;
% M = 10;
% Mx = 5;
% y = randn(nfft*M,2);
% x = randn(nfft*Mx,2);
% [abetaN,F,aGFTcrit] = aGFT(y,x,nfft,nfft,0.05);
% figure;plot(F,abetaN,'b',[0 F(end)],[aGFTcrit aGFTcrit],'k--')
% axis([0 F(end) 0 1]);xlabel('Frequency(Hz)');ylabel('MGFT')
% text(10,min(aGFTcrit*1.2,.98), ['FA = ' num2str(mean(double(abetaN>aGFTcrit))*100) '%'],'fontsize',12)
% fprintf('\nFA = %f\n', mean(double(abetaN>aGFTcrit)))
%---------------------------------------------------------------------

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
Svv = sum(Y(1:nf,:,:),2);
Snm = sum(X(1:nf,:,:),2);

abeta(1,:) = mean(Svv./(Svv+Snm),3);

if nargin > 3
    F = (0:(nf-1))*fs/tj; %frequency vector
    if nargin == 5
        nRuns = 10000;
        y  = betarnd(M,Mx,nRuns,N);
        aux = mean(y,2);
        abetacrit = quantile(aux,1-alpha);
    end
end

end