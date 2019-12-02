% CHMF Circular hybrid median filter.
%   J = CHMF(I,W) performs the circular hybrid median filter for despeckling
%   a breast ultrasound image I by using a window size WxW, where W is odd.
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   J = chmf(I,9);
%   figure; 
%   subplot 121; imshow(I); title('Original Image');
%   subplot 122; imshow(mat2gray(J)); title('Filtered Image');
%
%   See also ADF ADLG ISF ISFAD
%
%
%   Reference:
%   ----------
%   W. Gomez, W. C. A. Pereira, "A contrast enhancement method for 
%   improving the segmentation of breast lesions on ultrasonography," 
%   Computers in Medicine and Biology, vol.80, pp. 14-23, 2017.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   CHMF Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function f = chmf(g,m)
if mod(m,2) == 0
    error('Window size MxM must be odd');
end
% Creates circular windows for filtering
m2 = round(m/2);
[X,Y] = meshgrid(-m2+1:m2-1);
C = double(sqrt(X.^2 + Y.^2)<m2);
theta = zeros(m);
theta(X>0&Y>=0) = atan(Y(X>0&Y>=0)./X(X>0&Y>=0));
theta(X>0&Y<0) =  atan(Y(X>0&Y<0)./X(X>0&Y<0))+2*pi;
theta(X<0) = atan(Y(X<0)./X(X<0))+pi;
theta(X==0&Y>0) = pi/2;
theta(X==0&Y<0) = 3*pi/2;
theta = fix(theta*180/pi);
D = [22.5 67.5;112.5 157.5;202.5 247.5;292.5 337.5];
WD = zeros(m);
for i = 1:4
   WD(theta>=D(i,1) & theta<D(i,2)) = 1; 
end
WR = double(~WD); WR(m2,m2) = 0;
WR = WR.*C; WD = WD.*C;
% Filter image
f1 = ordfilt2(g, fix(median(1:sum(WD(:)))), WD, 'symmetric');
f2 = ordfilt2(g, fix(median(1:sum(WR(:)))), WR, 'symmetric');
M = cat(3,g,f1,f2);
f = fix(median(double(M),3));