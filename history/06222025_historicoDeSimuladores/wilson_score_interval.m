function [phat, pci] = wilson_score_interval(k, n, confidence_level)
% WILSON_SCORE_INTERVAL Calcula o intervalo de confiança de Wilson para uma proporção.
% k: número de sucessos
% n: número total de trials
% confidence_level: nível de confiança (ex: 0.95 para 95%)

if n == 0
    phat = NaN;
    pci = [NaN, NaN];
    return;
end

alpha = 1 - confidence_level;
z = norminv(1 - alpha/2); % Valor crítico z (ex: 1.96 para 95%)

phat = k/n;

% Componentes do intervalo de Wilson
term1 = phat + (z^2)/(2*n);
term2 = z * sqrt( (phat*(1-phat)/n) + (z^2)/(4*n^2) );
denominator = 1 + (z^2)/n;

pci_lower = (term1 - term2) / denominator;
pci_upper = (term1 + term2) / denominator;

% Garantir que os limites estejam entre 0 e 1
pci_lower = max(0, pci_lower);
pci_upper = min(1, pci_upper);

pci = [pci_lower, pci_upper];
end