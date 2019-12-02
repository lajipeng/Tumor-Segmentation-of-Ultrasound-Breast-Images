% CROPROI Crop a region of interest.
%   [J,K] = CROPROI(I,BW) crops a region of interest (ROI) of the gray-scale
%   image I by using the tumor limits of the binary image BW. The resultant
%   gray-scale ROI is in the image J and the corresponding binary shape in 
%   the image K.
%
%   [J,K] = CROPROI(I,BW,OPT) if OPT=1 forces the ROI to be a square matrix 
%   by cropping the largest dimension to make it equal to the lowest dimension.
%   By default OPT=0.
%
%   Example:
%   -------
%   load('BUS01.mat');
%   [J,K] = cropROI(I,Smanual);
%   delintumor(J,K);
%
%   See also DELINTUMOR

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   CROPROI Version 1.0 (Matlab R2014a Unix)
%   June 2017
%   Copyright (c) 2017, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [J,B] = cropROI(I,BW,opt)
if nargin < 3
    opt = 0;
end
[y,x] = find(BW);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
J = I(ymin:ymax,xmin:xmax);
B = BW(ymin:ymax,xmin:xmax);

if opt
    [M,N] = size(J);
    c = round(abs(M-N)/2);
    if N > M
       if ((mod(M,2)>0)&&(mod(N,2)>0))||((mod(M,2)==0)&&(mod(N,2)==0))
           J = J(:,c+1:N-c);
           B = B(:,c+1:N-c);
       elseif ((mod(M,2)>0)&&(mod(N,2)==0))||((mod(M,2)==0)&&(mod(N,2)>0))
           J = J(:,c:N-c);
           B = B(:,c:N-c);
       end
    elseif M > N
       if ((mod(M,2)>0)&&(mod(N,2)>0))||((mod(M,2)==0)&&(mod(N,2)==0))
           J = J(c+1:M-c,:);
           B = B(c+1:M-c,:);
       elseif ((mod(M,2)==0)&&(mod(N,2)>0))||((mod(M,2)>0)&&(mod(N,2)==0))
           J = J(c:M-c,:);
           B = B(c:M-c,:);
       end
    end
end