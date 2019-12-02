% PHASECONG Phase congruency of a gray-scale image.
%   PC = PHASECONG(I) calculates the measure of phase congruency of the 
%   gray-scale image I. By default, six scales and eight orientations are
%   used. PC is a cell array of phase congruency images for each orientation.
%
%   PC = PHASECONG(I,NSCALE,NORIENT) calculates the measure of phase congruency 
%   of the gray-scale image I with NSCALE scales and NORIENT orientations.
%
%   [PC,sPC] = PHASECONG(I,NSCALE,NORIENT) obtains the overall phase congruency 
%   results in PC array.
%
%   [PC,sPC,EO] = PHASECONG(I,NSCALE,NORIENT) obtains 2D cell array (NSCALExNORIENT)
%   of complex valued convolution results.  The real part is the result of convolving 
%   with the even symmetric filter, the imaginary part is the result from convolution 
%   with the odd symmetric filter.
%
%   NOTE: This function is based on the Peter Kovesi's implementation of phase  
%   congruency. School of Computer Science & Software Engineering. The University of 
%   Western Australia.
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   [PC,sPC] = phasecong(I);
%   figure; 
%   subplot 121; imshow(mat2gray(I));   title('Original');
%   subplot 122; imshow(mat2gray(sPC)); title('Overall PC');
%
%   See also MULTILOGGABOR MULTIRANKLET
%
%
%   References:
%   ----------
%   P. Kovesi, "Symmetry and asymmetry from local phase," In: Proceedings, 
%   10th Australian Joint Conference on Artificial Intelligence. 
%   Cambridge, MA: MIT Press; 1997, pp. 185-190.
%   
%   P. Kovesi, "Image Features From Phase Congruency," Videre: A Journal of 
%   Computer Vision Research, MIT Press, vol. 1, no. 3, 1999.
%
%   http://www.peterkovesi.com/matlabfns/PhaseCongruency/Docs/convexpl.html

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   PHASECONG Version 1.0 (Matlab R2014a Unix)
%   June 2017
%   Copyright (c) 2017, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [PC,pcSum,EO] = phasecong(im,nscale,norient)
% Parameters
if nargin < 2
    nscale          = 6;     % Number of wavelet scales.    
    norient         = 8;     % Number of filter orientations.
end
minWaveLength   = 3;     % Wavelength of smallest scale filter.    
mult            = 2;     % Scaling factor between successive filters.    
sigmaOnf        = 0.55;  % Ratio of the standard deviation
cutOff          = 0.4;   % Filter response spread
k               = 2.0;   % No of standard deviations of the noise
g               = 10;    % Gain factor of the sigmoid function
epsilon         = .0001;
% Fourier transform of image
im = double(im);
imagefft = fft2(im);	
% Image size
[rows,cols] = size(im);
zero = zeros(rows,cols);
EO = cell(nscale, norient);
PC = cell(norient,1);
EnergyV = zeros(rows,cols,3);
pcSum = zeros(rows,cols); 
% Set up X and Y matrices with ranges normalised to +/- 0.5
if mod(cols,2)
    xrange = (-(cols-1)/2:(cols-1)/2)/(cols-1);
else
    xrange = (-cols/2:(cols/2-1))/cols; 
end
if mod(rows,2)
    yrange = (-(rows-1)/2:(rows-1)/2)/(rows-1);
else
    yrange = (-rows/2:(rows/2-1))/rows; 
end
[x,y] = meshgrid(xrange, yrange);
% Normalized radius from center
radius = sqrt(x.^2 + y.^2); 
% Polar angle
theta = atan2(-y,x);        
% Quadrant shift radius and theta
radius = ifftshift(radius); 
theta  = ifftshift(theta);
radius(1,1) = 1;
sintheta = sin(theta);
costheta = cos(theta);
% Construct a low-pass filter that is as large as possible
lprad = sqrt(x.^2 + y.^2);
lp = ifftshift(1.0./(1.0 +(lprad./.45).^(2*15)));
logGabor = cell(1,nscale);
for s = 1:nscale
    wavelength = minWaveLength*mult^(s-1);
    fo = 1.0/wavelength;                  % Center frequency of filter
    logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));  
    logGabor{s} = logGabor{s}.*lp;        % Apply low-pass filter
    logGabor{s}(1,1) = 0;
end
for o = 1:norient 	% For each orientation...
    angl = (o-1)*pi/norient;
    ds = sintheta * cos(angl) - costheta * sin(angl);
    dc = costheta * cos(angl) + sintheta * sin(angl);
    dtheta = abs(atan2(ds,dc)); 
    dtheta = min(dtheta*norient/2,pi);
    spread = (cos(dtheta)+1)/2; 
    % Initialize accumulator matrices
    sumE_ThisOrient   = zero;          
    sumO_ThisOrient   = zero;       
    sumAn_ThisOrient  = zero;      
    Energy            = zero; 
    for s = 1:nscale
        lfilter = logGabor{s} .* spread;                   % Multiply radial and angular components
                                                           % to get a log-Gabor filter 
        EO{s,o} = ifft2(imagefft .* lfilter);              % Convolve image with even and odd filters
        An = abs(EO{s,o});                                 % Amplitude of even & odd filter response
        sumAn_ThisOrient = sumAn_ThisOrient + An;          % Sum of amplitude responses
        sumE_ThisOrient = sumE_ThisOrient + real(EO{s,o}); % Sum of even filter convolution results
        sumO_ThisOrient = sumO_ThisOrient + imag(EO{s,o}); % Sum of odd filter convolution results
        if s == 1 
            tau = median(sumAn_ThisOrient(:))/sqrt(log(4));    % Use median to estimate noise statistics
            maxAn = An;
        else
            maxAn = max(maxAn,An); % Maximum amplitude of components across scales
        end
    end
    % Energy vector data
    EnergyV(:,:,1) = EnergyV(:,:,1) + sumE_ThisOrient;
    EnergyV(:,:,2) = EnergyV(:,:,2) + cos(angl)*sumO_ThisOrient;
    EnergyV(:,:,3) = EnergyV(:,:,3) + sin(angl)*sumO_ThisOrient;
    % Get weighted mean filter response vector
    XEnergy = sqrt(sumE_ThisOrient.^2 + sumO_ThisOrient.^2) + epsilon;   
    MeanE = sumE_ThisOrient ./ XEnergy; 
    MeanO = sumO_ThisOrient ./ XEnergy; 
    % Calculate An(cos(phase_deviation) - | sin(phase_deviation)) | 
    for s = 1:nscale,       
        E = real(EO{s,o}); O = imag(EO{s,o});    % Extract even and odd convolution results.
        Energy = Energy + E.*MeanE + O.*MeanO - abs(E.*MeanO - O.*MeanE);
    end
    totalTau = tau * (1 - (1/mult)^nscale)/(1-(1/mult));
    EstNoiseEnergyMean = totalTau*sqrt(pi/2);        % Expected mean and std
    EstNoiseEnergySigma = totalTau*sqrt((4-pi)/2);   % values of noise energy
    T = EstNoiseEnergyMean + k*EstNoiseEnergySigma;  % Noise threshold
    Energy = max(Energy - T, 0);                     % Apply noise threshold
    % Calculate fractional 'width' of the frequencies
    width = (sumAn_ThisOrient./(maxAn + epsilon))/(nscale); 
    % Calculate the sigmoidal weighting function for this orientation.
    weight = 1.0 ./ (1 + exp(g*(cutOff - width)));
    % Apply weighting to energy and then calculate phase congruency
    PC{o} = weight.*Energy./sumAn_ThisOrient;
    % Overall phase congruency
    pcSum = pcSum+PC{o};
end