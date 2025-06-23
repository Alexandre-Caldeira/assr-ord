function res = funcoesPrimitivas(dna, sinal, sinalB, sinalLFT,~, M, fr)
%FUNCOES Executa uma função representada por uma árvore
%   Cada operador e função matemática utilizada é representado por um
%   índice que está contido na primeira posição do vetor data de cada nó da
%   árvore. Um fator de ponderamento ocupa a segunda posição deste vetor.
%   Esta função executa recursivamente cada nó a partir da raiz.
%
%   Entradas
%       sinal  - Matriz MxNxO das transformadas de fourier, apenas das
%                frequências positivas e já excluindo-se o nível DC. A
%                frequência de amostragem deve ser identica ao tamanho de
%                uma janela (para facilitar a simulação, já que o indice do
%                bin corresponderá à frequência do mesmo).
%                M - Bins
%                N - Janelas
%                O - Sujeitos
%       sinalLFT - Versão do sinal calculado sem janelamento.
%       sinalB - Matriz identica ao sinal, porem sempre contera apenas
%                ruido, equivale a uma coleta sem estimulação
%       M      - Número de janelas
%       fr     - Frequência a ser testada

L = 10;         % Número de bins laterais utilizados no teste F local

    % Função que executa a função descrita pela arvore combinatoria contida no DNA
        if dna==0                  %  0 - Constante
            res = 1;
        elseif dna==1              %  1 - MSC
            res = msc();
        elseif dna==2              %  2 - CSM
            res = csm();
        elseif dna==3              %  3 - TFG
            res = tfg();
        elseif dna==4              %  4 - TFL
            res = tfl();
        elseif dna==5              %  5 - MSC numerador 
            res = msc_n();
        elseif dna==6              %  6 - MSC denominador
            res = msc_d();
        elseif dna==7              %  7 - CSM parcela seno
            res = csm_s();
        elseif dna==8              %  8 - CSM parcela cosseno
            res = csm_c();
        elseif dna==9              %  9 - TFG numerador
            res = tfg_n();
        elseif dna==10             % 10 - TFG denominador
            res = tfg_d();
        elseif dna==11             % 11 - TFL numerador
            res = tfl_n();
        elseif dna==12             % 12 - TFL denominador
            res = tfl_d();
        elseif dna==13             % 13 - CSM hiperbólica
            res = csm_h();
        elseif dna==14             % 14 - TFL média
            res = tfl_media();
        end
        
    %% MSC
    % Versão completa
    function r = msc()
        r = (abs(sum(sinal(fr,1:M,:))).^2)./(sum(abs(sinal(fr,1:M,:)).^2))./M;
    end
    % Numerador
    function r = msc_n()
        r = (abs(sum(sinal(fr,1:M,:))).^2)./M;
    end
    % Denominador
    function r = msc_d()
        r = (sum(abs(sinal(fr,1:M,:)).^2));
    end

    %% CSM
    % Versão completa
    function r = csm()
        r = (sum(cos(angle(sinal(fr,1:M,:))))./M).^2+(sum(sin(angle(sinal(fr,1:M,:))))./M).^2;
    end
    % Parcela do seno
    function r = csm_s()
        r = (sum(sin(angle(sinal(fr,1:M,:))))./M).^2;
    end
    % Parcela do cosseno
    function r = csm_c()
        r = (sum(cos(angle(sinal(fr,1:M,:))))./M).^2;
    end
    % Versão hiperbólica
    function r = csm_h()
        r = abs(1-var(cos(angle(sinal(fr,1:M,:))),[],2)./mean(std(cos(angle(sinal(fr+[-5:-1 1:5],1:M,:))),[],2),1));
    end
    
    %% Teste F local normalizado
    % Versão completa
    function r = tfl()
        r = (abs(mean(sinalLFT(fr*M,1,:))).^(2))./sum(abs(mean(sinalLFT((fr*M)+[(-L:-L/2)+4 0 (L/2:L)-4],1,:),2)).^(2));
    end
    % Numerador
    function r = tfl_n()
        r = (abs(mean(sinalLFT(fr*M,1,:))).^2);
    end
    % Denominador
    function r = tfl_d()
        r = sum(abs(mean(sinalLFT((fr*M)+[-L:-L/2 L/2:L],1,:),2)).^2)./L;
    end
    % Versão adaptada onde o denominador inclui a média dos bins vizinhos e não a soma pura
    function r = tfl_media()
        r = (abs(mean(sinalLFT(fr*M,1,:))).^2)./(sum(abs(mean(sinalLFT((fr*M)+[-L/2:-1 0 1:L/2],1,:),2)).^2)./L+(abs(mean(sinalLFT(fr*M,1,:))).^2));
    end
    
    %% Teste F Global normalizado
    % Versão completa
    function r = tfg()
        r = sum(abs(sinal(fr,1:M,:)).^2)./...
            (sum(abs(sinal(fr,1:M,:)).^2)+sum(abs(sinalB(fr,1:M,:)).^2));
    end
    % Numerador
    function r = tfg_n()
        r = sum(abs(sinal(fr,1:M,:)).^2);
    end
    % Denominador
    function r = tfg_d()
        r = (sum(abs(sinal(fr,1:M,:)).^2)+sum(abs(sinalB(fr,1:M,:)).^2));
    end
end