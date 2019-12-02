% HISTEQU Histogram equalization.
%   J = HISTEQU(I) performs the histogram equalization to enhance the 
%   contrast of a 8-bit image I.
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   J = histequ(I);
%   figure; 
%   subplot 121; imshow(I); title('Original Image');
%   subplot 122; imshow(J); title('Contrasted Image');
%
%   See also CLAHE FUZZYENH SACE SIGMOIDFILT  
%
%
%   Reference:
%   ----------
%   R.C. Gonzalez, R.E. Woods, Digital Image Processing, Pearson Prentice Hall, 
%   New Jersey, USA, 2008.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   HISTEQU Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function J = histequ(I)
I = double(I);
N = numel(I);
% Image histogram
Hc = accumarray(I(:)+1,ones(N,1),[256 1],@sum,0);
% Cumulative sum
c = 255/N;
fa = c.*cumsum(Hc);
% Transformation function
J = uint8(fa(I+1));