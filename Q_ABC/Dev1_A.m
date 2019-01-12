clear;
close all;
clc;
nvec = [10,100,1000];
m = 5000;

% Algorithme
for n = nvec;
    
    Xvec = [];
    for i = 1:m
        X = 0;
        for j = 1:n
            Saut = round(rand(1));
            if Saut == 0;
                Saut = -1;
            end
            X = X+Saut;
        end
        Xvec = [Xvec,X];
    end
    
    % Moyenne et Écart-type
    Moyenne = mean(Xvec);
    Std = std(Xvec);
    
    % Histogramme
    figure;
    hist(Xvec,(max(Xvec)-min(Xvec))/2+1)
    xlabel('X')
    ylabel('Occurence')
    title(strcat('n = ',num2str(n)))
    txt = sprintf(strcat('Moyenne = ',num2str(Moyenne),'\nSTD = ',num2str(Std)));
    TextLocation(txt,'Location','NW')
end