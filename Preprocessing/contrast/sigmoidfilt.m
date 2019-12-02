% SIGMOIDFILT Sigmoid filter.
%   J = SIGMOIDFILT(I,ALPHA,BETA) performs the sigmoid filter to enhance the 
%   contrast of a breast ultrasound image I, where ALPHA is the width
%   parameter and BETA the center parameter.
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   J = sigmoidfilt(I,8,60);
%   figure; 
%   subplot 121; imshow(I); title('Original Image');
%   subplot 122; imshow(J); title('Contrasted Image');
%
%   See also CLAHE FUZZYENH HISTEQU SACE  
%
%
%   Reference:
%   ----------
%   W.K. Moon, S.-C. Chang, C.-S. Huang, R.-F. Chang, "Breast tumor classification 
%   using fuzzy clustering for breast elastography," Ultrasound in Medicine 
%   and Biology, vol. 37, no. 5, pp. 700-708, 2011.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   SIGMOIDFILT Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------


function J = sigmoidfilt(I,a,b)
I = double(I);
Imax = max(I(:));
Imin = min(I(:));
J = uint8(((Imax-Imin)./(1+exp(-(I-b)/a)))+Imin);