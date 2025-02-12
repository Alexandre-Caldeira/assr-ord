%Função que converte o valor do detector para o seu respectivo p_value na
%condição de HO.
%Entrada: 
% s=> string (nome do detector) 
% x=> valor do detector
% parametrs => parametros da distribuição na condição Ho.
%   parameters.M => número de janelas 
%   parameters.L =>  banda lateral (teste-F)
%Saída
%  pvalue=> p value 

%Exemplo:
%1) MSC
%M=50; alfa = 0.05; kcrit = 1-alfa^(1/(M-1));
%parameters.M = M;
%pvalue = det2pvalue('MSC',kcrit,parameters) 

% CSM
%M=50; alfa = 0.05;kcrit = chi2inv(1-alfa,2)/(2*M);
%parameters.M = M;
%pvalue = det2pvalue('CSM',kcrit,parameters) 
%LFT
%L=24; alfa = 0.05;kcrit = finv((1-alfa),2,2*L);
%parameters.L = L;
%pvalue = det2pvalue('LFT',kcrit,parameters) 

function pvalue = det2pvalue(s,x,parameters) 



%% Detector
if strcmp(s,'MSC') == 1 
   b = (parameters.M-1);
   pvalue = betacdf(x,1,b);
 %  pvalue = betapdf(x,1,b);
end

if strcmp(s,'CSM') == 1 
   x= x*2*parameters.M;
   pvalue = chi2cdf(x,2);
end

if strcmp(s,'LFT') == 1 
   pvalue = fcdf(x,2,2*parameters.L);
end


