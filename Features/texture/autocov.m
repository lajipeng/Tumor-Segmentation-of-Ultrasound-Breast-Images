% AUTOCOV Gray-level autocovariance coefficients.
%   [X,FEAT] = AUTOCOV(I,BW) computes 24 the autocovariance coefficients 
%   from the gray-level image I masked by the binary image BW. X is a numeric 
%   vector with the feature values and FEATS is a cell vector with the name of 
%   the features in the same order as in X.
%
%   Example:
%   -------
%   load('BUS01.mat');   
%   [x,feats] = autocov(I,Smanual);
%
%   See also AUTOCORR BDIP_BVLC GLCM LAWSENERGY
%
%   Reference:
%   ---------
%   R.-F. Chang, W.-J. Wu, W.K. Moon, D.-R. Chen, "Improvement in breast 
%   tumor discrimination by support vector machines and speckle-emphasis 
%   texture analysis," Ultrasound in Medicine and Biology, vol. 29, no. 5,
%   pp. 679-686 2003.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   AUTOCOV Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = autocov(US,BW)
[y,x] = find(BW);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
I = double(US(ymin:ymax,xmin:xmax));
[M,N] = size(I);
f_bar = mean2(I); % Media imagen
% Computa matriz de autocovarianza de 5x5
sz = 5;
A = zeros(sz);
fea = cell(sz);
for dM = 1:sz
   for dN = 1:sz
        i = 1:(M - dM);
        j = 1:(N - dN);
        F1 = I(i,j)-f_bar;
        F2 = I(i+dM,j+dN) - f_bar;
        A(dM,dN) = (1/((M-dM)*(N-dN)))*sum(sum(F1.*F2));
        fea(dM,dN) = {['ACOV_(' num2str(dN-1) ',' num2str(dM-1) ')']};
   end
end
% Salidas
N = (A./A(1,1))';
N = N(:)'; N(1) = [];
fea = fea'; fea = fea(:); fea(1) = [];
feats = fea';
x = N;