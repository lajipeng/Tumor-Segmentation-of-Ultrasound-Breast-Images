% PAB Posterior acoustic behavior.
%   [X,FEATS] = PAB(I,BW) computes two posterior acoustic behavior features
%   from the gray-level image I masked by the binary image BW: posterior acoustic 
%   feature class (PS_D) and minimum side difference (MSD). X is a numeric vector 
%   with the feature values and FEATS is a cell vector with the name of the features 
%   in the same order as in X.
%
%   [X,FEATS] = PAB(BW,FEATURES) computes a specific set of texture
%   features according with FEATURES. The set of valid strings includes 
%   (case insensitive):
%
%       'all'   - Entire feature set.
%       'psd'   - Posterior acoustic feature class.
%       'msd'   - Minimum side difference.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS01.mat');   
%   [x,feats] = pab(I,Smanual);
%
%   Example 2: Compute the minimum side difference
%   ----------------------------------------------
%   load('BUS01.mat');   
%   [x,feats] = pab(I,Smanual,'msd');
%
%   See also AVMASS NRG
%
%
%   References:
%   ----------
%   W.-C. Shen, R.-F. Chang, W.K. Moon, Y.-H. Chou, C.-S. Huang, "Breast ultrasound
%   computer-aided diagnosis using bi-rads features," Academic Radioliology,
%   vol. 14, no. 8, pp. 928-939, 2007.
%
%   K. Horsch, M. L. Giger, L. A. Venta, C. J. Vyborny, "Computerized diagnosis 
%   of breast lesions on ultrasound," Medical Physics, vol. 29, no. 2, 
%   pp. 157-164, 2002.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   PAB Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = pab(US,BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'psd';'msd'}; 
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
US = double(US);
BW = logical(BW);
[N,M] = size(US);
[y,x] = find(BW);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
gap = round((xmax-xmin)/6);%esto para quitar un sexto de cada lado
% Posterior Acoustic Feature Class
UScrop = US(ymin:N,xmin+gap:xmax-gap);
BWcrop = BW(ymin:N,xmin+gap:xmax-gap);
BWl = bwlabel(~BWcrop);
[N2,M2] = size(UScrop);
BWpos = BWl==BWl(N2,round(M2/2));
d = sum(BWpos);
if any(d>100)
   a = max(d)-100;
   BWpos = BWpos(1:N2-a,1:M2);
   UScrop = UScrop(1:N2-a,1:M2);
end
PS_d = mean(UScrop(BWpos)) - mean(US(BW));
% Minimum side difference
USpos = US(ymax:N,xmin+gap:xmax-gap);
USlft = US(ymax:N,1:xmin-1);
USrgt = US(ymax:N,xmax+1:M);
MSD = min(median(USpos(:))-median(USlft(:)),median(USpos(:))-median(USrgt(:)));
% salidas
aux_x = [PS_d,MSD];
aux_f = {'PS_D','MSD'};
x = aux_x(idxStats);
feats = aux_f(idxStats);