function X = espgram(x,tj)

%Fourier Transform of a signal divided in windows 
%
%X = espgram(signal,tj) 
%
%Parameters:  
%
%signal => Vector containing the signal for analysis
%tj => number of points of each segment




[l,c] = size(x); %Transform x[n] in column vector
if l ==1
    x = x(:);
end

[tamsinal,c] = size(x);
ncol = fix(tamsinal/tj); %number of windows 
x = x(1:ncol*tj,1); %limit the signal for a entire number of windows. Prevents shorter windows.
x = reshape(x,tj,ncol); %transform the signal in a matrix where each column is a retangular window of the signal  

if 0
for i=1:ncol,
    x(:,i) = x(:,i).*blackman(tj); %to be used when other type of windowing is required, e.g. for hanning window, use hann(tj) 
end
end

X = fft(x); %Fast Fourier Transform

nfft = fix(tj/2); %number of points in the returned Fourier Transform
X = X(1:nfft+1,:); %only half of the FFT-values are returned, because of FFT's simmetry in real signals

% Antonio Mauricio - somewhere in time
% Leo Bonato - 13/11/2003
% Leo Bonato - 17/11/2003 => Wrote the help
