function output = getShotNoise(p, photonCounts,...
    pixelWidth, widthPSFs, bkgdSigma, plotFLG)
%
%
%
%
%
% Edited:
%   KGryte - (2012-05-14) - Created.
%
% References:
%   Holden et al (2011) Biophys J.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Checks!

if nargin < 6
    
    plotFLG = false; % default;
    
end % end IF;

if nargin < 3
    
    % Pixel Width:
    pixelWidth = 94 * 10^-9; % nm; 'a'; measured on our EMCCD camera

    % Standard deviation of background photon counts per pixel:
    bkgdSigma = sqrt([2.9, 2.9]); % photons; 'b_{d}', 'b_{a}'; [1 x numChannels]
    
    % Widths of channel PSFs:
    widthPSFs = [132*10^-9. 150*10^-9]; % nm; 's_{d}', 's_{a}'; [1 x numChannels]
   
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization:

% Excess noise factor:
noiseFactor = sqrt(2); % 'f_{g}'






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shot Noise Calculation:

% Compute the first quantity of our shot noise calculation: (this is the
% 'classic' shot noise equation for FRET, with the addition of an excess noise factor
% from the EMCCD camera)
X1 = noiseFactor^2 .* p .* (1-p) / photonCounts;

% Calculate the channel photon counts, using the ratiometric equation for
% E:
c2 = round(p .* photonCounts); % 'A'; 0 < x < 1
c1 = round((1-p) .* photonCounts); % 'D'

% Compute the second quantity of the shot noise calculation: (see Holden et
% al (2011); his equation is wrong, need to switch c1 and c2 in his eqn.)

scalar = (4*pi)/(pixelWidth^2 .* photonCounts.^4);

term1 = c2.^2 * bkgdSigma(1)^2 * widthPSFs(1)^2;
term2 = c1.^2 * bkgdSigma(2)^2 * widthPSFs(2)^2;

X2 = scalar .* (term1 + term2); 

% Now compute the shot noise:
shotNoise = sqrt(X1 + X2); % For low photon counts, this behaves as a smiley, rather than a frowny face. noiseFactor^2*



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot:

if plotFLG
    
    plot(x, shotNoise, 'k');
    
end % end IF



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output:


output = shotNoise;



                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EOF

