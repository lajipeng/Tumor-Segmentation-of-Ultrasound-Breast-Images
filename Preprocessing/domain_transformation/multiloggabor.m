% MULTILOGGABOR Multichannel decomposition by Log-Gabor filters.
%   J = MULTILOGGABOR(I,NSCALES,NORIENTATIONS,OCTAVES,Q) performs the multichannel 
%   decomposition of input image I by a bank of Log-Gabor filters, where NSCALES
%   is an integer to indicate the number of scales, NORIENTATIONS is an integer
%   to indicate number of orientations, Q is the quantization level of
%   output images, and OCTAVES is the bandwidth of the filters in octaves
%   and could take one of the following values: 0.5, 1, 1.5 and 2. All the 
%   orientations at the same scale are averaged; hence, if the  size of the input 
%   image I is NxM then the output J is a hypermatrix of size MxNxNSCALES.
%
%   [J,LG] = MULTILOGGABOR(I,NSCALES,NORIENTATIONS,OCTAVES,Q) obtains the
%   individual log-Gabor channels in the cell array LG, where the columns 
%   represent the scales and the rows represent the orientations. No quantization 
%   is applied to the channels.
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   J = multiloggabor(I,6,12,2,32);
%   figure; 
%   subplot 231; imshow(mat2gray(J(:,:,1))); title('Scale 1');
%   subplot 232; imshow(mat2gray(J(:,:,2))); title('Scale 2');
%   subplot 233; imshow(mat2gray(J(:,:,3))); title('Scale 3');
%   subplot 234; imshow(mat2gray(J(:,:,4))); title('Scale 4');
%   subplot 235; imshow(mat2gray(J(:,:,5))); title('Scale 5');
%   subplot 236; imshow(mat2gray(J(:,:,6))); title('Scale 6');
%
%   See also MULTIRANKLET PHASECONG QUANTIZATION
%
%
%   References:
%   ----------
%   P. Kovesi, "Symmetry and asymmetry from local phase," In: Proceedings, 
%   10th Australian Joint Conference on Artificial Intelligence. 
%   Cambridge, MA: MIT Press; 1997, pp. 185-190.
%
%   W. Gomez, W. C. A. Pereira, A. F. C. Infantosi, "Breast ultrasound 
%   despeckling using anisotropic diffusion guided by texture descriptors," 
%   Ultrasound in Medicine and Biology, vol. 40, no. 11, pp. 2609-2621, 2014.
%
%   W. Gomez, B. A. Ruiz-Ortega, "New fully-automated method for segmentation
%   of breast lesions on ultrasound based on texture analysis," Ultrasound
%   in Medicine and Biology, vol. 42, no. 7, pp. 1637-1650, 2016.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   MULTILOGGABOR Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [PGab,GC] = multiloggabor(I,nscale,norient,octave,q)

% sigmaOnf  .90   mult 1.15
% sigmaOnf  .85   mult 1.2
% sigmaOnf  .75   mult 1.4       (bandwidth ~1 octave)
% sigmaOnf  .65   mult 1.7
% sigmaOnf  .55   mult 2.2       (bandwidth ~2 octaves)

if octave == 0.5
    sigmaOnf = .85;
    mult = 1.2;
elseif octave == 1
    sigmaOnf = .75;
    mult = 1.4;
elseif octave == 1.5
    sigmaOnf = .65;
    mult = 1.7;    
elseif octave == 2
    sigmaOnf = .55;
    mult = 2.2;
else
    error('Invalid octave!');
end

[Mi,Ni] = size(I);
LGab = LogGaborFilt(I, nscale, norient, 3, mult, sigmaOnf);
PGab = zeros(Mi,Ni,nscale);
GC = cell(norient,nscale);
for i = 1:nscale
    aux = zeros(Mi,Ni,norient);
    for j = 1:norient
        lg = LGab{i,j};
        aux(:,:,j) = lg;
        GC{j,i} = lg;
    end
    if nargin == 5
        PGab(:,:,i) = quantization(mean(aux,3),q);
    else
        PGab(:,:,i) = mean(aux,3);
    end
end

%*************************************************************************
%*************************************************************************
% Kovesi's implementation of Log-Gabor filters
% Copyright (c) 2001-2010 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.peterkovesi.com/matlabfns/PhaseCongruency/Docs/convexpl.html
function LG = LogGaborFilt(IO, nscale, norient, minWaveLength, mult, sigmaOnf)

dThetaOnSigma = 1.3;
TF = fft2(IO); 
TF(1,1) = 0;

[rows, cols] = size(IO);
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
radius = sqrt(x.^2 + y.^2);
theta = atan2(y,x);
radius = ifftshift(radius);
theta  = ifftshift(theta);
radius(1,1) = 1;
sintheta = sin(theta);
costheta = cos(theta);
clear x; clear y; clear theta;
lp = lowpassfilter(xrange,yrange,.4,10);
logGabor = cell(1,nscale);
for s = 1:nscale
    wavelength = minWaveLength*mult^(s-1);
    fo = 1.0/wavelength;
    logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));  
    logGabor{s} = logGabor{s}.*lp;
    logGabor{s}(1,1) = 0;
    L = sqrt(sum(logGabor{s}(:).^2));
    logGabor{s} = logGabor{s}./L;        
end
% The main loop
LG = cell(nscale,norient);
for o = 1:norient
    angl = (o-1)*pi/norient;
    wavelength = minWaveLength;
    ds = sintheta * cos(angl) - costheta * sin(angl);
    dc = costheta * cos(angl) + sintheta * sin(angl);
    dtheta = abs(atan2(ds,dc));
%     dtheta = min(dtheta*norient/2,pi);
%     spread = (cos(dtheta)+1)/2;        
    thetaSigma = pi/norient/dThetaOnSigma;  
    spread = exp((-dtheta.^2) / (2 * thetaSigma^2));  
    for s = 1:nscale
        filter = logGabor{s} .* spread;
        L = sqrt(sum(real(filter(:)).^2 + imag(filter(:)).^2 ))/sqrt(2);
        LG(s,o) = {real(ifft2(TF.*(filter./L)))};
        wavelength = wavelength * mult;
    end
end
%------------------------------------------------------------------------
function f = lowpassfilter(xrange, yrange, cutoff, n)
[x,y] = meshgrid(xrange, yrange);
radius = sqrt(x.^2 + y.^2);
f = ifftshift( 1.0 ./ (1.0 + (radius ./ cutoff).^(2*n)) );
%*************************************************************************
%*************************************************************************