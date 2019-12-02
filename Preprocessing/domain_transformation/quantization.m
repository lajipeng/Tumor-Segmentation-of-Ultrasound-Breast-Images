% QUANTIZATION Gray-level quantization.
%   J = QUANTIZATION(I,NL,GL) performs the quantization of the input image I
%   to NL gray levels. GL is a a two-element vector, [LOW HIGH], that 
%   specifies how the values in I are scaled into gray levels. If GL is not 
%   given or is empty QUANTIZATION uses the minimum and maximum grayscale 
%   values in I as limits, [min(I(:)) max(I(:))].
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   J = quantization(I,8);
%   figure; 
%   subplot 121; imshow(I); title('Original Image');
%   subplot 122; imshow(mat2gray(J)); title('Quantized to 8 levels');
%
%   See also MULTILOGGABOR MULTIRANKLET
%
%
%   Reference:
%   ---------
%   W. Gomez, W.C.A. Pereira, A.F.C. Infantosi, "Analysis of co-occurrence 
%   texture statistics as a function of gray-level quantization for classifying 
%   breast ultrasound," IEEE Trans. Med. Imaging, vol. 31, no. 10, 
%   pp. 1889-1899, 2012.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   QUANTIZATION Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function SI = quantization(I,NL,GL)
I = double(I);
if (nargin < 3)||isempty(GL)
    GL = [min(I(:)) max(I(:))];
end
slope = NL / (GL(2) - GL(1));
intercept = 1 - (slope*(GL(1)));
SI = round(imlincomb(slope,I,intercept,'double'));
SI(SI > NL) = NL;
SI(SI < 1) = 1;