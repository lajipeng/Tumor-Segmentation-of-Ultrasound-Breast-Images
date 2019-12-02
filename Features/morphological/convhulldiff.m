% CONVHULLDIFF Morphological features based on convex hull.
%   [X,FEATS] = CONVHULLDIFF(BW) computes three morphological features from the
%   binary blob of a breast lesion BW: overlap ratio, normalized residual value,
%   and convexity. X is a numeric vector with the feature values and FEATS 
%   is a cell vector with the name of the features in the same order as in X.
%
%   [X,FEATS] = CONVHULLDIFF(BW,FEATURES) computes a specific set of 
%   morphological features according with FEATURES. The set of valid 
%   strings includes (case insensitive):
%
%       'all'   - Entire feature set.
%       'or'    - Overlap ratio.
%       'nrv'   - Normalized residual value.
%       'cnvx'  - Convexity.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = convhulldiff(BW);
%   % Same as [x,feats] = convhulldiff(BW,'all')
%
%   Example 2: Compute the overlap ratio
%   ------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = convhulldiff(BW,'or');
%
%   Example 3: Compute convexity and normalized residual value
%   ----------------------------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = convhulldiff(BW,'cnvx','nrv');
%
%   See also EQUIVELLIPSE GEOMETRIC NRL
%
%   
%   References:
%   ----------
%   A.V. Alvarenga, A.F.C. Infantosi, W.C.A. Pereira, C.M. Azevedo, "Assessing 
%   the performance of morphological parameters in distinguishing breast tumors 
%   on ultrasound images," Med. Eng. Phys., vol. 32, no. 1, pp. 49-56, 2009.
%
%   R.-F. Chang, W.-J. Wu, W. Moon, D.-R. Chen, "Automatic ultrasound segmentation
%   and morphology based diagnosis of solid breast tumors," Breast Cancer Res. Treat.,
%   vol. 89, no. 2, pp. 179-185, 2005.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   CONVHULLDIFF Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = convhulldiff(BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'OR';'NRV';'CNVX'};      
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
BW = double(BW);
BW_props = regionprops(BW,'ConvexHull','Perimeter','Centroid');
CH = roipoly(BW,BW_props.ConvexHull(:,1),BW_props.ConvexHull(:,2));
OR  = bwarea(BW)/bwarea(CH);         % Overlap ratio
NRV = bwarea(xor(BW,CH))/bwarea(CH); % Normalized Residual Value
CH_props = regionprops(CH,'Perimeter');
Cnv = CH_props.Perimeter/BW_props.Perimeter; % Convexity
aux_x = [OR NRV Cnv];
aux_f = {'OveRatio','NRV','Convexity'};
% Outputs
x = aux_x(idxStats);
feats = aux_f(idxStats);
