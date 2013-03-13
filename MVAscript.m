%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script for analyzing RPo data.
%       All species. 
%
%
%
%
%   Edited:
%       2012-07-12  -   KGryte: Created.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear mvaobj
clear results

results{8} = [];

% Filename:
fileName = 'Particle'; 

% Get the particle numbers
numPart = [1,2,3,4,5,6,8,...
    10,11,12,13,15,16,17,18,19,...
    20,21,23,26,29,...
    33,35,36,37,38,39,...
    40,41,42,43,44,46,47,48,...
    54,55,57,59,...
    65,...
    72,77,78,79,...
    81,85,89,...
    93,94,95,96,...
    102,104,105,106,108];

% Effective Gain:
U = 6;

pixelWidth = 94*10^-9;
widthPSFs = [132, 150] .* 10^-9; 
bkgdMean = [2.9, 2.9];

% Initialize the average counts vector:
avgCounts = nan(numel(numPart), 1);

dat = [];

% Loop through each particle and extract the relevant data:
for i = 1 : numel(numPart)
    
    % Create the dynamic file name:
    file = [fileName, int2str(numPart(i)), '.mat'];
    
    % Load the file:
    load(file, 'PlotData');
    
    % Get the relevant data:
    DexDem = PlotData.UnfilteredDexDem(1:100,2);
    DexAem = PlotData.UnfilteredDexAem(1:100,2);
    
    % Get the average counts:
    avgCounts(i) = round(nanmean(DexDem + DexAem) / U);
    
    % Calculate FRET:
    E = DexAem ./ (DexDem + DexAem);
    
    % Assign over the FRET to the variable dat:
    dat = [dat; E];
        
    % Perform MVA:
    mvaobj{i} = mva(E);
    mvaobj{i} = mvaobj{i}.analyzeData();
    
    for j = 1 : numel(mvaobj{i}.windowSizes)
        % Concatenate the results for each window size:
        results{j} = [results{j}; mvaobj{i}.results{j,1}];
        
    end % end FOR
      
    if mod(i,5)
        fprintf('%d ', i);
    else
        fprintf('%d \n', i);
    end % end IF/ELSE
    
end % end FOR i

% Calculate the mean photons:
meanPhotons = round(mean(avgCounts));

% Initialize a figure:
hFig = figure;

% Configure the figure colormap:
load('mva_colormap.mat', 'mva_colormap');

set(hFig,...
    'Colormap', mva_colormap);            

% Determine the number of subplots:
numSubplots = numel(mvaobj{1}.windowSizes);

% Create the edges:
xEdges = 0 : 1/100 : 1;
yEdges = logspace(-4, -1, 100+1); % linspace(0, .01, obj.nBins+1);      

offset = 1.30^2; % account for OLS fitting error
edgesSN = 0:0.01:1;

% Generate the histograms:
for i = 1 : numSubplots

    % Create the axes:
    hAx(i,1) = subplot(3,ceil(numSubplots/3),i);

    % Histogram the data:
    [Counts, edges, binCenters] = histcn(...
        results{i}(:,1:2),...
        xEdges,...
        yEdges); % Default is 100 bins

    % Smooth the histogram counts:
    smoothedData = smoothn(...
        Counts,...
        20);%smoothing); 

    % Create the contour plot:
    [ContourMatrix, hContour] = contourf(...
        binCenters{1},... % bin centers
        binCenters{2},... % bin centers
        smoothedData(:,:,1)',... Counts',...
        10,...
        'EdgeColor','none');

    set(hAx(i,1), 'YScale', 'log');                

    % Plot guidelines for the population mean and variance: 
    hold on 
    
%     globalMean = nanmean(results{i}(:,1));
%     globalVar = nanstd(results{i}(:,1)).^2;
% 
%     plot([0, 1], repmat(globalVar, 1, 2), 'k:'); % population variance
%     plot(repmat(globalMean, 1, 2), [min(yEdges), max(yEdges)], 'k:'); % population mean 
    
    shotNoise = getShotNoise(edgesSN, meanPhotons, pixelWidth, widthPSFs, sqrt(bkgdMean),  0); 
    
    plot(hAx(i), edgesSN, offset*shotNoise.^2, 'k'); %% Additional factor!

    hold off

    % Set the title:
    title(hAx(i,1), ['Window Size: ', int2str(mvaobj{1}.windowSizes(i))]);

    % Set the labels:
    xlabel(hAx(i,1), 'Mean');
    ylabel(hAx(i,1), 'Variance');


end % end FOR





% Make a histogram:
edges = 0:0.01:1;

counts = histc(dat, edges);

figure;
bar(edges, counts, 'histc');

xlim([0,1]);

hold on

shotNoiseHist = zeros(numel(edges),1);

for i = 1 : numel(numPart)
    
    % Build the shot noise distribution:
    photonCounts = ones(numel(mvaobj{i}.data), 1) .* avgCounts(i);
    
    % Include background in shot noise distribution generation:
    shotNoiseHist = shotNoiseHist + generateShotNoiseDistribution(photonCounts, nanmean(mvaobj{i}.data), edges, bkgdMean, pixelWidth, widthPSFs);
end 

stairs(edges, shotNoiseHist, 'r');

xlabel('Efficiency');
ylabel('Counts');

legend('Data', 'Shot-Noise');
legend boxoff

hold off



save('Analysis - MVA - (2012-07-12)(1).mat', 'mvaobj', 'results', 'dat');











return;

numData = 300;

U = 6; % Gain factor
DexDem = round(PlotData.UnfilteredDexDem(1:numData,2) ./ U);
DexAem = round(PlotData.UnfilteredDexAem(1:numData,2) ./ U);
figure; 
subplot(2,1,1); 
plot(...
    PlotData.UnfilteredDexDem(1:numData,1), DexDem, 'g',...
    PlotData.UnfilteredDexAem(1:numData,1), DexAem, 'r'); 
ylabel('Inferred Photons');
xlim([0,30])
subplot(2,1,2); 
plot(...
    PlotData.UnfilteredDexDem(1:numData,1),...
    PlotData.UnfilteredDexAem(1:numData,2) ./ (PlotData.UnfilteredDexDem(1:numData,2) + PlotData.UnfilteredDexAem(1:numData,2)), 'r');

xlabel('Time [sec]');
ylabel('Efficiency');
ylim([0,1]);
xlim([0,30])









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOS