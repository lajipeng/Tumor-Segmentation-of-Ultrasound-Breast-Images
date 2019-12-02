% MARGCLASS Margin class features.
%   [X,FEATS] = MARGCLASS(BW) computes two morphological features from the
%   binary blob of a breast lesion BW: undulation characteristic and angular
%   characteristic. X is a numeric vector with the feature values and FEATS 
%   is a cell vector with the name of the features in the same order as in X.
%
%   [X,FEATS] = MARGCLASS(BW,FEATURES) computes computes a specific set 
%   of morphological features according with FEATURES. The set of valid 
%   strings includes (case insensitive):
%
%       'all'   - Entire feature set.
%       'und'   - Undulation characteristic.
%       'ang'   - Angular characteristic.
%       'u+a'   - Undulation plus angular characteristics.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = margclass(BW);
%   % Same as [x,feats] = margclass(BW,'all');
%
%   Example 2: Compute angular characteristic
%   -----------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = margclass(BW,'ang');
%
%   See also NSPD_LI SPICULATION
%   
%
%   Reference:
%   ---------
%   W.-C. Shen, R.-F. Chang, W.K. Moon, Y.-H. Chou, C.-S. Huang, "Breast ultrasound
%   computer-aided diagnosis using bi-rads features," Academic Radioliology,
%   vol. 14, no. 8, pp. 928-939, 2007.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   MARGCLASS Version 2.0 (Matlab R2014a Unix)
%   December 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = margclass(BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'und';'ang';'u+a'};      
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
% Crop
[y,x] = find(BW);
xmn = min(x);
xmx = max(x);
ymn = min(y);
ymx = max(y);
BW2 = BW(ymn:ymx,xmn:xmx);
BW2 = padarray(BW2,[1 1],0,'both');
% Compute distance map
D = bwdist(~BW2);
% Inscribed circle distance map
[M,N] = size(BW2);
[y,x] = find(D==max(D(:)));
xc = mean(x); yc = mean(y);
[X,Y] = meshgrid(1:N,1:M);
C = sqrt((X-xc).^2 + (Y-yc).^2);
% Get lobules with the maximum inscribed circle
r = max(D(:)); % Radius
BW3 = xor(~BW2,C>=r);
D2 = bwdist(~BW3);
BW3(D2<2) = 0;
% Count undulations without considerinh slighter undulations
L = bwlabel(BW3);
nl = max(L(:));
U = 0;
for i = 1:nl
    idx = L==i;
    D3 = bwdist(~idx); % Interior distances
    if max(D3(:)) > 3
        U = U+1; % Its a lobule
    end
end
% Compute angularity feature
% Skeleton
S = bwmorph(BW2,'skel','inf'); % skeleton(BW2, D, 10, 1); %
BW4 = bwareaopen(and(S,~(C<r)),5);
L2 = bwlabel(BW4);
A = max(L2(:)); % Angularity
MUA = U+A;
% Features
aux_x = [U A MUA];
aux_f = {'mMU','mMA','mMUA'};
x = aux_x(idxStats);
feats = aux_f(idxStats);