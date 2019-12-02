% NRG Normalized radial gradient.
%   [X,FEAT] = NRG(I,BW) computes the normalized radial gradient feature from
%   the gray-level image I masked by the binary image BW. X is the numeric value 
%   of the feature and FEAT is the name of the feature.
%
%   Example:
%   -------
%   load('BUS01.mat');   
%   [x,feat] = nrg(I,Smanual);
%   disp(['x = [' num2str(x) ']']);
%   disp(feat);
%
%   See also AVMASS PAB
%
%
%   Reference:
%   ---------
%   K. Horsch, M. L. Giger, L. A. Venta, C. J. Vyborny, "Computerized diagnosis 
%   of breast lesions on ultrasound," Medical Physics, vol. 29, no. 2, 
%   pp. 157-164, 2002.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   NRG Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feat] = nrg(US,BW)
I = double(US);
[M,N] = size(I);
% 5x5 sobel filter
sx = [-1 -2  0 2  1
      -4 -8  0 8  4
      -6 -12 0 12 6
      -4 -8  0 8  4
      -1 -2  0 2  1];
sy = sx';
gx = imfilter(I,sx,'replicate');
gy = imfilter(I,sy,'replicate');
[y,x] = find(BW);
xc = mean(x);
yc = mean(y);
[x,y] = meshgrid(1:N,1:M);
x = x - xc;
y = y - yc;
n = sqrt(x.^2+y.^2);
xu = x./n;
yu = y./n;
F1 = gx.*xu + gy.*yu;
F2 = sqrt(gx.^2+gy.^2);
C = bwperim(BW);
x = sum(F1(C))/sum(F2(C));
feat{1} = 'nrg';