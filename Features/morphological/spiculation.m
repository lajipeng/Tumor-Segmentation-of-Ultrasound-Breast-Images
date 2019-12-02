% SPICULATION Spiculation features of the lesion.
%   [X,FEATS] = SPICULATION(BW) computes two spiculation features from the
%   binary blob of a breast lesion BW: spiculation and number of skeleton
%   end-points. X is a numeric vector with the feature values and FEATS 
%   is a cell vector with the name of the features in the same order as in X.
%
%   [X,FEATS] = SPICULATION(BW,FEATURES) computes a specific set of 
%   morphological features according with FEATURES. The set of valid 
%   strings includes (case insensitive):
%
%       'all'   - Entire feature set.
%       'spic' 	- Spiculation feature.
%       'skend' - Number of skeleton end-points.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = spiculation(BW);
%   % Same as [x,feats] = spiculation(BW,'all')
%
%   Example 2: Compute the number of skeleton end-points
%   ----------------------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = spiculation(BW,'skend');
%
%   See also NSPD_LI MARGCLASS
%
%
%   Reference:
%   ---------
%   Y. Su, Y. Wang, J. Jiao, Y. Guo, "Automatic detection and classification 
%   of breast tumors in ultrasonic images using texture and morphological 
%   features," Open Med. Inf. J., vol. 5, pp. 26-37, 2011.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   SPICULATION Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feat] = spiculation(BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'spic';'skend'};
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
BW = double(BW);
Prps = regionprops(BW,'Centroid');
xc = Prps.Centroid(1);
yc = Prps.Centroid(2);
Baux = bwboundaries(BW);
C = Baux{1};
C(:,1) = C(:,1)-yc;
C(:,2) = C(:,2)-xc;
[~,r] = cart2pol(C(:,2),C(:,1));
R = fft(r); N = length(R);
w = 2*pi*(0:(N-1))/N;
spic = sum(abs(R(w<=pi/4)))/sum(abs(R(pi/4<w&w<=pi/2)));
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
skl = bwmorph(BW,'skel',Inf);
skend = bwmorph(skl,'endpoints');
skl_end = sum(skend(:));
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
aux_x = [spic skl_end];
aux_f = {'Spicul' 'SkelEnds'};
x = aux_x(idxStats);
feat = aux_f(idxStats);