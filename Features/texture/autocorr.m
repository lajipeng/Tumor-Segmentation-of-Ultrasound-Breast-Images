% AUTOCORR Gray-level autocorrelation coefficient.
%   [X,FEAT] = AUTOCORR(I,BW) computes the autocorrelation feature from the
%   gray-level image I masked by the binary image BW. X is the numeric value 
%   of the feature and FEAT is the name of the feature.
%
%   Example:
%   -------
%   load('BUS01.mat');   
%   [x,feat] = autocorr(I,Smanual);
%
%   See also AUTOCOV BDIP_BVLC GLCM LAWSENERGY   
%
%
%   Reference:
%   ---------
%   K. Horsch, M. L. Giger, L. A. Venta, C. J. Vyborny, "Computerized diagnosis 
%   of breast lesions on ultrasound," Medical Physics, vol. 29, no. 2, 
%   pp. 157-164, 2002.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   AUTOCORR Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feat] = autocorr(US,BW)
% Rectangulo minimo
[y,x] = find(BW);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
K = double(US(ymin:ymax,xmin:xmax));
R = (K/max(K(:))).^2;
[NR,MR] = size(R);
% Autocorrelacion en profundidad
Cy = zeros(NR-1,MR);
for n = 1:NR-1
	a = R(n+1:NR,:);
	b = R(1:NR-n,:);
	Cy(n,:) = sum(a.*b,1);
end
% Autocorrelacion lateral
Cy_bar = sum(Cy,2);
% Autocorrelacion total
x = sum(Cy_bar./Cy_bar(1));
feat = {'AutoCorr'};