clc; clear; clf; close all;

trialCount = 3;
file = "1500f";



for trial = 1:1:trialCount
    name = file + int2str(trial) + "_data.txt";
    data = readmatrix(name);

    if trial == 1
        fullData = data;
        N = length(data(:,1));
    else
        if(length(data(:,1)) ~= N)
            nShort = length(data(:,1));
            inN = N - nShort;
            if(inN < 0)
                fullData = [fullData; data(N+1:end, 1), NaN(abs(inN), trial-1)];
            else
                data(:,2) = [data(:,2); NaN(inN, 1)];
            end
        end
        fullData = [fullData, data(:, 2)];
    end 
end 


f = figure(1);
f.Theme = "Light";

hold on 

t = fullData(:,1);

for trial = 1:1:trialCount
    dataName = "Trial " + int2str(trial);
    plot(t, fullData(:,trial + 1), '-o', 'DisplayName', dataName)
end 



meanData = mean(fullData(1:N, 2:trialCount + 1), 2);
plot(t(1:N), meanData, '-o', "LineWidth", 4, "DisplayName", "Mean of All Trials")

legend()
xlabel('Time (s)');
ylabel('Angular Displacement (deg? - confirm units)');
title('Trial Data Comparison');

