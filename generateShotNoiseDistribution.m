function shotNoiseHist = generateShotNoiseDistribution(numTrials, ...
    probability, histEdges, bkgdMean, pixelWidth, widthPSFs)
% GENERATESHOTNOISEDISTRIBUTION     Generates the shot noise distribution
% for the ratio of a binomial random variable and the total number of
% trials.
%
%   Shot noise describes the statistical fluctations resulting from stochastic
%   processes.
%
%
%   Parameters:
%       numTrials:      a vector of trial numbers (e.g., photon counts)
%
%       probability:    the success probability for a binomial random
%                       variable (e.g., transfer efficiency between a donor
%                       and acceptor fluorophore; FRET)
%       histEdges:      edges defining the bins into which the shot noise
%                       predicted values should fall (e.g., 0:0.01:1)
%
%   Output:
%       shotNoiseHist:  vector of histogram counts describing the
%                       distribution of successes (in terms of the proportion
%                       of number of trials; e.g., how many sets of trials
%                       had a success percentage between 40 and 50%?)
%
%
%
%   Dependencies:
%       binornd.m
%       randn.m
%
%   Edited:
%       2012-07-05: KGryte - Created.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization:

% Oversampling factor: (addresses infrequently sampled trial numbers)
oversamplingFactor = 100;

% Data length:
numData = size(numTrials(:),1); % number of rows in the column vector

% Initialize the shot noise distribution array:
shotNoiseDistr = nan(numData, oversamplingFactor);


% Determine if background parameters provided:
if nargin > 3
    includeBkgd = true;
else
    includeBkgd = false;
end % end IF/ELSE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shot Noise Distribution:


if ~includeBkgd

    % Run through each photon count:
    for i = 1 : numData

        % [1] For N trials, we draw a random variable A from a binomial
        % distribution with parameters N and p:
        randVals = binornd(numTrials(i), probability, 1, oversamplingFactor); % [1 x oversamplingFactor] vector

        % [2] Calculate the ratio between A and N:
        ratio = randVals ./ numTrials(i);

        % [3] Add the ratio values to our shot noise distribution:
        shotNoiseDistr(i, :) = ratio; % insert a row vector into the array 

    end % end FOR


    
else % include background when generating shot noise distribution
    
    % Run through each photon count:
    for i = 1 : numData
        
        % [1] Determine the background counts:
        bkgdTerm = 4.*pi.*(widthPSFs./pixelWidth).^2 .* bkgdMean; % [1 x numChannels]; Term 2 from Equation 19, Thompson et al
        
        %
        for m = 1 : oversamplingFactor
            
            % Generate new background for each sampling:
            totalBkgd = round(bkgdTerm + (sqrt(bkgdTerm).*randn)); % assume st dev = sqrt(mean); distribution approximated by Gaussian; [1 x numChannels]

            % No negative backgrounds:
            totalBkgd(totalBkgd < 0) = 0;

            % [2] For N trials, we draw a random variable A from a binomial
            % distribution with parameters N and p:
            N = numTrials(i) - sum(totalBkgd); % by reducing N (remove photons contributed by background), we increase the dispersion of the ratio.
            randVal = binornd(N, probability, 1, 1); %oversamplingFactor); % [1 x oversamplingFactor] vector

            % [3] Calculate the ratio between A and N:
            ratio = randVal ./ N;

            % [4] Add the ratio value to our shot noise distribution:
            shotNoiseDistr(i, m) = ratio; % insert into a row vector into the array 
            
        end % end FOR m

    end % end FOR
    
    
end % end IF/ELSE




% Compute the histogram:
counts = histc(shotNoiseDistr(:), histEdges);

% Divide the counts by the oversampling factor:
shotNoiseHist = counts ./ oversamplingFactor;
    
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOF