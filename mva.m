classdef mva
    % the mean-variance analysis class
    %
    %
    %
    %
    %
    %   Dependencies:
    %       - BUFFER.M: Signal Processing toolbox
    %
    %
    %   Edited:
    %       KGryte - (2012-05-14) - Created.
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % ------------------------------------------------------------------- %
    %% PROPERTIES:
    
    properties(GetAccess = 'public', SetAccess = 'public')
        % Public properties:
        windowSizes = [3,5,10,15,20,30,50,100]; % in data points
        data;
        results;
        globalMean;
        globalVar; 
        
    end % end PROPERTIES
    
    
    properties(GetAccess = 'protected', SetAccess = 'protected')
        % Private properties:
        nBins = 100;
        smoothing = 20;
        xMin = 0;
        xMax = 1;
        yMin = -4; % 10e-4
        yMax = -1; % 10e-1    
        xEdges;
        yEdges;
        
    end % end PROPERTIES
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        % Restricted properties:
        hFig;
        hAx;
       
    end % end PROPERTIES
    
    % ------------------------------------------------------------------- %
    %% METHODS:
    
    
    % ------------------------ %
    % Public Methods:
    methods(Access = 'public')
        
        %% CONSTRUCTOR:
        function obj = mva(data, windowSizes)
            % Class constructor:
            
            if nargin > 0
                
                obj.data = data;
                
                if nargin > 1
                    obj.windowSizes = windowSizes;
                end % end IF
                
            end % end IF
            
                
                
            
        end % end function MVA
        
        %% METHOD: analyzeData
        function obj = analyzeData(obj)
            % ANALYZEDATA 
            
            for i = 1 : numel(obj.windowSizes)
                
                % Grab the window size:
                W = obj.windowSizes(i);
                
                % Use a function from the Signal Processing toolbox: buffer.m
                %   - this allows us to create an array of overlapping
                %   windows from our data set:
                dat = buffer(obj.data, W, W-1, 'nodelay'); % [WxN-W+1]
                
                % Calculate the sample mean and sample variance
                meanVals = obj.sampleMean(dat); % row vector [1xN-W+1]
                varVals = obj.sampleVariance(dat,repmat(meanVals, W, 1)); % row vector [1xN-W+1]
                
                % Assign over to RESULTS:
                obj.results{i, 1} = [meanVals',  varVals']; % [N-W+1x2]
                
                
            end % end FOR    
            
            % Calculate the global mean and variance: (global here refers
            % to the population mean and variance, which are effectively
            % larger samples of the parent distribution)
            obj.globalMean = obj.sampleMean(obj.data);
            obj.globalVar = obj.sampleVariance(obj.data, obj.globalMean);
            
            
        end % end function ANALYZEDATA
        
        %% METHOD: plotData
        function obj = plotData(obj)
            % PLOTDATA
            
            % Initialize a figure:
            obj.hFig = figure;
            
            % Configure the figure colormap:
            load('mva_colormap.mat', 'mva_colormap');

            set(obj.hFig,...
                'Colormap', mva_colormap);            
            
            % Determine the number of subplots:
            numSubplots = numel(obj.windowSizes);
            
            % Create the edges:
            obj.xEdges = obj.xMin : (obj.xMax-obj.xMin)/obj.nBins : obj.xMax;
            obj.yEdges = logspace(obj.yMin, obj.yMax, obj.nBins+1); % linspace(0, .01, obj.nBins+1);             
            
            % Generate the histograms:
            for i = 1 : numSubplots
                
                % Create the axes:
                obj.hAx(i,1) = subplot(3,ceil(numSubplots/3),i);
                
                % Histogram the data:
                [Counts, edges, binCenters] = histcn(...
                    obj.results{i,1}(:,1:2),...
                    obj.xEdges,...
                    obj.yEdges); % Default is 100 bins
                
                % Smooth the histogram counts:
                smoothedData = smoothn(...
                    Counts,...
                    obj.smoothing); 

                % Create the contour plot:
                [ContourMatrix, hContour] = contourf(...
                    binCenters{1},... % bin centers
                    binCenters{2},... % bin centers
                    smoothedData(:,:,1)',... Counts',...
                    10,...
                    'EdgeColor','none');
                                
                set(obj.hAx(i,1), 'YScale', 'log');                
                
                % Plot guidelines for the population mean and variance: 
                hold on 
                
                plot([obj.xMin, obj.xMax], repmat(obj.globalVar, 1, 2), 'k:'); % population variance
                plot(repmat(obj.globalMean, 1, 2), [min(obj.yEdges), max(obj.yEdges)], 'k:'); % population mean 
                
                hold off
                
                % Set the title:
                title(obj.hAx(i,1), ['Window Size: ', int2str(obj.windowSizes(i))]);
                
                % Set the labels:
                xlabel(obj.hAx(i,1), 'Mean');
                ylabel(obj.hAx(i,1), 'Variance');
                
                
            end % end FOR
            
        end % end function PLOTDATA
            
        
    end % end METHODS (PUBLIC)
    
    
    % ----------------------- %
    % Private Methods:
    methods(Access = 'private');
        
        %% METHOD: sampleMean
        function val = sampleMean(obj, data)
            % SAMPLEMEAN calculate the sample mean
            
            % Calculate the sum:
            sumData = sum(data,1); % row vector
            
            % Divide by the sample size:
            val = sumData ./ size(data,1); % row vector
            
        end % end function SAMPLEMEAN
        
        %% METHOD: sampleVariance
        function val = sampleVariance(obj, data, mu)
            % SAMPLEVARIANCE calculate the sample variance
            
            % Subtract the mean value from each data point:
            val = data - mu;
            
            % Square each value:
            val = val.^2;
            
            % Sum the values:
            val = sum(val,1); % row vector
            
            % Divide by 1 less than the sample size:
            val = val ./ (size(data,1) - 1); % row vector
            
        end % end function SAMPLEVARIANCE
        
        
    end % end METHODS (PRIVATE)
    
    
end % end CLASSDEF



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOF