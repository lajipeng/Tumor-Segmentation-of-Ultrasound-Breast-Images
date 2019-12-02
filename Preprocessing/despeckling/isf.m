% ISF Interference-based speckle filter.
%   J = ISF(I,W) performs the interference-based speckle filter for 
%   despeckling a breast ultrasound image I by using a window size WxW, 
%   where W is odd. 
%   
%   Note: The ISF filter is based on the circular hybrid median filter
%   instead of the circular median filter originally proposed  by 
%   Cardoso et al.
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   J = isf(I,9);
%   figure; 
%   subplot 121; imshow(I); title('Original Image');
%   subplot 122; imshow(mat2gray(J)); title('Filtered Image');
%
%   See also ADF ADLG CHMF ISFAD
%
%
%   Reference:
%   ----------
%   F. M. Cardoso, M. M. S. Matsumoto, S. S. Furuie "Edge-preserving speckle
%   texture removal by interference-based speckle filtering followed by 
%   anisotropic diffusion," Ultrasound in Medicine and Biology, vol. 38,
%   no. 8, pp. 1414-1428, 2012.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   ISF Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function J = isf(I,w)
I = double(I);
Imed = chmf(I,w);	% Circular hybrid median filter
Ic   = max(I,Imed);	% Destructive interference suppression
J    = chmf(Ic,w);  % Circular hybrid median filter