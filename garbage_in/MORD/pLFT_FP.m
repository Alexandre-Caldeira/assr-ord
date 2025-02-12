function [FyN, FyNcrit] = pLFT_FP(y,L,fs,alpha,fo)
%
% Additive local spectral F-test.
% 
% Extension of the Local F-test to multi channel analysis as the sum of 
% local spectral F-test of each channel. This is a MORD technique that
% makes use only magnitude information and can be used as a detector of
% hidden periodicities in noise, since N real-valued signals are available.
% The method is described in "Felix, L.B., Rocha, P.F.F., Mendes, E.M.A.M.
% and Miranda de Sá, A.M.F.L.. Multivariate approach for estimating the
% local spectral F-test and its application to the EEG during photic
% stimulation. Computer Methods and Programs in Biomedicine, p. 87-91,
% 2018."
%
% Sintaxes:
% [FyN, FyNcrit] = ALFT(y,fo,L,fs) returns estimated value of Additive
% Local F-test at fo and respective simulated critical value.
%
% Inputs:
% y => matrix whose columns are the signals to be considered. 
% fo => frequency to be tested
% L => number of closest neighbouring frequencies to fo. L must be even.
% fs => sample rate.
%
% Example
%
% fs=100;M=10;fo=21;L=38;
% y1 = awgn(sin(2*pi*fo*(linspace(0,M,M*fs))),-15,'measured','db')';
% y2 = awgn(sin(2*pi*fo*(linspace(0,M,M*fs))),-16,'measured','db')';
% Rfo = fo-(L/2):1/M:fo+(L/2);
% MLphi=zeros(size(Rfo));
% for i=1:length(Rfo)
% [MLphi(i),MLphicrit] = aLFT([y1 y2],L,fs,0.05,Rfo(i));
% end
% figure;plot(Rfo,MLphi,'b',fo,1.05*MLphi(fix(length(Rfo)/2)+1),'-rv',[Rfo(1) Rfo(end)],[MLphicrit MLphicrit],'k--')
L=12; 
if bitget(L,1) %odd
    clc
    L = L+1;
    disp('L must be even. Now L=L+1');
end

[t,N] = size(y);   % length and number of channel
nfft = fix(t/2);   % Number of points for FFT
pfo = round(fo*nfft/(fs/2))+1; % Position FFT bin of fo

% Amplitude spectrum in fo and closest neighbouring frequencies
Y = fft(y);
Y = abs(Y(1:nfft+1,:)); 
Yfo = Y(pfo,:); 
Yfn(1:L,1:N) = [Y(pfo-1,:);Y(pfo-2,:);Y(pfo-3,:);Y(pfo-4,:);Y(pfo-8,:);Y(pfo-9,:);Y(pfo+1,:);Y(pfo+2,:);Y(pfo+3,:);Y(pfo+4,:);Y(pfo+8,:);Y(pfo+9,:)]; 

% Yfn = Y(pfo-L/2:pfo+L/2,:);
% Yfn(L/2+1,:) = []; 

% Compute F value
FyN = (prod((Yfo.^2)./(1/L*sum(Yfn.^2,1)))).^(1/N);

% % Critical Value
% load('SimFyNcrit.mat');
% FyNcrit = Fcrit(N,L/2); % for alpha=0.05

if nargout>1
    [FyNcrit]=pLFTcrit(N,L,alpha);
end

end

function [pLFTcrit]=pLFTcrit(N,L,alpha)
% aLFT critical values using Monte Carlo
%
% 
nRuns = 1e6; %number of iterations. It runs instantaneuosly in a i5-7600 CPU @ 3.5Hz (mar/2019)
y =  frnd(2,2*L,nRuns,N);
pLFT = prod(y,2); 
pLFTcrit = quantile(pLFT,1-alpha);
end

