% AVMASS Gray-level average-based texture features.
%   [X,FEATS] = AVMASS(I,BW) computes three average-based texture features
%   from the gray-level image I masked by the binary image BW: lesion boundary
%   class (abrupt interface), average gray intensity, and echo pattern feature. 
%   X is a numeric vector with the feature values and FEATS is a cell vector 
%   with the name of the features in the same order as in X.
%
%   [X,FEATS] = AVMASS(BW,FEATURES) computes a specific set of texture
%   features according with FEATURES. The set of valid strings includes 
%   (case insensitive):
%
%       'all'   - Entire feature set.
%       'lb'    - Lesion boundary.
%       'epi'   - Average internal gray intensity.
%       'epg'   - Echo pattern feature.
%       'epc'   - Contrast feature.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS01.mat');   
%   [x,feats] = avmass(I,Smanual);
%   % Same as [x,feats] = avmass(I,Smanual,'all')
%
%   Example 2: Compute lesion boundary feature
%   ------------------------------------------
%   load('BUS01.mat');   
%   [x,feats] = avmass(I,Smanual,'lb');
%
%   Example 3: Compute lesion boundary and echo pattern features
%   ------------------------------------------------------------
%   load('BUS01.mat');   
%   [x,feats] = avmass(I,Smanual,'lb','epg');   
%   
%   See also NRG PAB
%
%   
%   Reference:
%   ---------
%   W.-C. Shen, R.-F. Chang, W.K. Moon, Y.-H. Chou, C.-S. Huang, "Breast ultrasound
%   computer-aided diagnosis using bi-rads features," Academic Radiology,
%   vol. 14, no. 8, pp. 928-939, 2007.
%
%    W.-C. Shen, R.-F. Chang, W.K. Moon, "Computer aided classification system for 
%    breast ultrasound based on breast imaging reporting and data system 
%   (bi-rads)", Ultrasound in Medicine and Biology, vol. 33, no. 11, 
%   pp. 1688-1698, 2007.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   AVMASS Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = avmass(US,BW,varargin)
%*********************************************************************
% Check for options
if nargin < 3
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'LB';'EPi';'EPg';'EPc'};      
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
US = double(US);
BW = logical(BW);
% Abrupt interface
k = 3;
Din  = bwdist(~BW);
Dout = bwdist(BW);
Win  = Din<=k&Din>0;
Wout = Dout<=k&Dout>0;
LB = mean(US(Wout))-mean(US(Win));
% Average gray intensity
EPi = mean(US(BW));
% Internal echo pattern feature
hx = [1 2 1;0 0 0;-1 -2 -1];
hy = hx';
Gx = imfilter(US,hx);
Gy = imfilter(US,hy);
G = sqrt(Gx.^2+Gy.^2);
EPg = mean(G(BW));
% Internal pattern
J = US(BW)+1;
N = sum(BW(:));
H = accumarray(J,ones(N,1),[256 1],@sum,0);
p = 0.25*N;
for i = 256:-1:1
    if sum(H(i:256)) > p
        break;
    end
end
EPc = (mean(J(J>=i)) - EPi)/EPi;
% Salidas
aux_x = [LB,EPi,EPg,EPc];
aux_f = {'LB_d','EP_i','EP_g','EP_c'};
x = aux_x(idxStats);
feats = aux_f(idxStats);