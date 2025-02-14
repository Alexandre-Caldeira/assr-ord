function [aThresholds, gThresholds] = vc_beta_cgst(M, K)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

TotalAlpha      = 0.05;                        
Alpha_k         = ones(1,K)*(TotalAlpha/K);     
Gamma_k         = ((1-TotalAlpha)/K).*ones(1,K);

Resolution      = (1/0.0001);                  
Xvalues         = 0:1/Resolution:3;            
Null         	= betapdf(Xvalues, 1, M-1); %normpdf(Xvalues,muk1,sigmak1);% normpdf(Xvalues,0,0.05); %normpdf(Xvalues,0.05,0.3); %betapdf(Xvalues, 1, M-1); 
Null            = Null/sum(Null);             	
Beta_Norm       = Null/sum(Null);             	
k               = 1;                            
aThresholds(k)	= 1 - Alpha_k(k).^(1./(M-1));  % quantile(betarandn(1,1e5), 1-Alpha_k(1)); %
gThresholds(k)	= 1-(1- Gamma_k(k)).^(1./(M-1)); % quantile(randn(1,1e5), Gamma_k(1));
TruncInd_Ra     = round(aThresholds(k)*Resolution);                                         
TruncInd_Rg     = round(gThresholds(k)*Resolution);           

for k = 2:K
    % disp(k)
    NullTrunc                   = Null;                                                     
    NullTrunc(TruncInd_Ra:end)  = zeros(1, length(NullTrunc(TruncInd_Ra:end)));              
    NullTrunc(1:TruncInd_Rg)    = zeros(1, length(NullTrunc(1:TruncInd_Rg)));
    Null2                       = conv(Beta_Norm, NullTrunc);                              
    Null2                       = Null2 / (sum(Null2) / (1 - sum(Gamma_k(1:(k-1))) - sum(Alpha_k(1:(k-1)))));
    TruncInd_Ra                 = findIndex(Null2, sum(Null2) - Alpha_k(k));            
    aThresholds(k)              = TruncInd_Ra/Resolution;                                    
    TruncInd_Rg                 = findIndex(Null2, Gamma_k(k), 1);
    gThresholds(k)              = TruncInd_Rg/Resolution;
    Null                        = Null2;
end

end