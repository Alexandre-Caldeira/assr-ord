function [ORD] = msc_fft(Y_input, M_win)
% MSC_FFT Calcula a Magnitude Squared Coherence (MSC).
% Y_input é a FFT dos dados, com dimensão [Nfreqs, M_windows]
% M_win é o número de janelas (M_windows)

% (O restante do código da função msc_fft conforme definido anteriormente)
% ... por exemplo ...

if size(Y_input,2) == 0 % Sem janelas, sem MSC
    ORD = zeros(size(Y_input,1),1);
    return;
end
if M_win == 0
    ORD = zeros(size(Y_input,1),1);
    return;
end

current_M_win = M_win;
if current_M_win > size(Y_input,2)
     current_M_win = size(Y_input,2);
     if current_M_win == 0
         ORD = zeros(size(Y_input,1),1);
         return;
     end
end

if current_M_win == 1 % Caso de uma única janela, MSC é 1.
    ORD = ones(size(Y_input,1),1);
    return;
end

numerator = abs(sum(Y_input(:, 1:current_M_win), 2)).^2;
denominator = current_M_win * sum(abs(Y_input(:, 1:current_M_win)).^2, 2);

ORD = numerator ./ denominator;
ORD(denominator == 0) = 0; % Se não há energia, MSC é 0
ORD(isnan(ORD)) = 0;
ORD(isinf(ORD)) = 0;
end % Certifique-se de que a função termina com 'end'