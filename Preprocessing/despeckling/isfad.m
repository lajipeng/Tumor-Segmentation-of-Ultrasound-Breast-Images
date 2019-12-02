% ISFAD Interference-based speckle filter followed by anisotropic diffusion.
%   J = ISFAD(I,W,MAXITER) performs the interference-based speckle filter 
%   followed by anisotropic diffusion for despeckling a breast ultrasound 
%   image I, where W is the window size of the median filter, and 
%   MAXITER is the maximum number of iterations for diffusion.
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   J = isfad(I,10,500);
%   figure; 
%   subplot 121; imshow(I); title('Original Image');
%   subplot 122; imshow(mat2gray(J)); title('Filtered Image');
%
%   See also ADF ADLG CHMF ISF
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
%   ISFAD Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function J = isfad(I,r,niter)
% Validations
if nargin < 2
    r = 10;
    niter = 500;
end
if nargin < 3
    niter = 500;
end
I = double(I);
Imed = medfilt(I,r); % Circular median filter
Ic = max(I,Imed);    % Destructive interference suppresion
ISF = medfilt(Ic,3); % Circular median filter
J = adf(ISF,1,0.5,niter,2); % Anisotropic diffusion filter
%**********************************************************************
function J = medfilt(I,r)
% Circular window
[x,y] = meshgrid(-r:r);
h = double(sqrt(x.^2+y.^2)<=r);
% Filter image
J = ordfilt2(I, fix(median(1:sum(h(:)))), h, 'symmetric');