function [r,F,rcrit] = csm(y,tj,fs,P)

%[r,F,rcrit] = csm(y,tj,fs,P);
%
%Component Synchrony Measure of a signal. Can be used as a detector of
%periodic signal in noise.
%
%CSM varies from 0 to 1 and estimates the degree to which the phases of the 
%frequency of interest are dispersed (CSM = 0 means phases are evenly dispersed
%across a 360° range among the subaverages) or clustered (PC = 1 means phases 
%are identical in all subaverages). 
% 
%Sintaxes:
%
%[r] = csm(y,tj) returns only the csm of the signal
%[r,F] = csm(y,tj,fs) returns CSM and frequency vector 
%[r,F,rcrit] = csm(y,tj,fs,P) for all aforementioned output plus the 
%associated critical value
%
%Input Parameters:
%
%y => signal whose phase clustering will be analized  
%tj => number of points of each segment for FFT routine
%fs => sample rate
%P => significance level, e.g. P = 0.05 

[l,~] = size(y); %Transform y[n] in column vector
if l ==1
    y = y(:);
end


[Y] = espgram(y,tj); %Fourier Transform of a signal divided in windows 

[~,M] = size(Y);
t = angle(Y); %phase angle 

r=(1/M^2)*sum(cos(t')).^2+(1/M^2)*sum(sin(t')).^2; %CSM of the signal

if nargin > 2
    F = 0:fs/tj:fs/2; %frequency vector
    if nargin == 4
        rcrit = chi2inv(1-P,2)/(2*M); % critical value
    end
end

%Antonio Mauricio - somewhere in time
%Leo Bonato - 18/11/2003