% GEOMETRIC Geometric features of the lesion.
%   [X,FEATS] = GEOMETRIC(BW) computes four morphological features from the
%   binary blob of a breast lesion BW: form factor, roundness, extent, and
%   depth-to-width ratio. X is a numeric vector with the feature values
%   and FEATS is a cell vector with the name of the features in the same
%   order as in X.
%
%   [X,FEATS] = GEOMETRIC(BW,FEATURES) computes a specific set of 
%   morphological features according with FEATURES. The set of valid 
%   strings includes (case insensitive):
%
%       'all'   - Entire feature set.
%       'form' 	- Form factor.
%       'round' - Roundness.
%       'ext'   - Extent.
%       'dwr'   - Depth-to-width ratio.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = geometric(BW);
%   % Same as [x,feats] = geometric(BW,'all')
%
%   Example 2: Compute the roundness
%   --------------------------------
%   load('BUS02.mat');   
%   [x,feats] = geometric(BW,'round');
%
%   Example 3: Compute extent and form factor
%   -----------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = geometric(BW,'ext','form');
%
%   See also CONVHULLDIFF EQUIVELLIPSE NRL
%
%
%   References:
%   ----------
%   C.-M. Chen, Y.-H. Chou, K.-C. Han, et al., "Breast lesions on sonograms: 
%   computer-aided diagnosis with nearly setting-independent features and 
%   artificial neural networks," Radiology, vol. 226, no. 2, pp. 504-514,2003.
%
%   R.-F. Chang, W.-J. Wu, W. Moon, D.-R. Chen, "Automatic ultrasound segmentation 
%   and morphology based diagnosis of solid breast tumors," Breast Cancer Res. Treat.,
%   vol. 89, no. 2, pp. 179-185, 2005.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   GEOMETRIC Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = geometric(BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'form';'round';'ext';'dwr'};      
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
BW = double(BW);
S = regionprops(double(BW),'Area','Perimeter','MajorAxisLength');
peri_BW = S.Perimeter;
area_BW = S.Area;
max_diameter = S.MajorAxisLength;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FormFactor = (4*pi*area_BW)/(peri_BW^2);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Roundness  = (4*area_BW)/(pi*max_diameter^2);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Ext = extent(BW);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[yBW,xBW] = find(BW);
xBWmax = max(xBW); xBWmin = min(xBW);
yBWmax = max(yBW); yBWmin = min(yBW);
DWR = (yBWmax-yBWmin)/(xBWmax-xBWmin);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
aux_x = [FormFactor Roundness Ext DWR];
aux_f = {'FormFactor' 'Roundness' 'Extent' 'D:W'};
x = aux_x(idxStats);
feats = aux_f(idxStats);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function ext = extent(BW)
A = bwarea(BW);
[x,y] = find(bwperim(BW));
x = x-mean(x);
y = y-mean(y);
a = 90;
t = (0:a-1)*pi/180;
c = cos(t);
s = sin(t);
xr = x*c - y*s;
yr = x*s + y*c;
Dh = abs(max(xr,[],1)-min(xr,[],1));
Dv = abs(max(yr,[],1)-min(yr,[],1));
ext = max(A./(Dh.*Dv));