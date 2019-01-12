clear;
close all;
clc;
nvec = [10,100,1000];
m = 5000;

% Algorithme
for n = nvec;
    
    Rvec = [];
    for i = 1:m
        X = 0;
        Y = 0;
        Z = 0;
        for j = 1:n
            Saut = 2*rand(1)-1;
            X = X+Saut;
        end
        for j = 1:n
            Saut = 2*rand(1)-1;
            Y = Y+Saut;
        end
        for j = 1:n
            Saut = 2*rand(1)-1;
            Z = Z+Saut;
        end
        Rvec = [Rvec,sqrt(X^2+Y^2+Z^2)];
    end
    
    % Moyenne et Écart-type
    Moyenne = mean(Rvec);
    Std = std(Rvec);
    Std_norm = std(Rvec)/(sqrt(n));
    
    % Histogramme
    figure;
    hist(Rvec,(max(Rvec)-min(Rvec))*2+1)
    xlabel('R')
    ylabel('Occurence')
    title(strcat('n = ',num2str(n)))
    txt = sprintf(strcat('Moyenne = ',num2str(Moyenne),'\nSTD = ',num2str(Std),'\nSTD/(n^{0.5}) = ',num2str(Std_norm)));
    TextLocation(txt,'Location','NE')
end