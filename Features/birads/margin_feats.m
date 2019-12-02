% MARGIN_FEATS Compute BI-RADS margin features.
%   [X,FEAT] = MARGIN_FEATS(BW) computes margin BI-RADS features according 
%   to the BI-RADS lexicon defined for breast ultrasound, where BW is the 
%   binary shape of the lesion:
%   
%   BI-RADS feature         Quantitative feature
%   ---------------         ----------------------------  
%   Margin
%                           Number of undulations (U)
%                           Angularity feature (A)
%                           Sum of U and A
%   
%   Example:
%   -------
%   load('BUS01.mat');   
%   [x,feat] = margin_feats(Smanual);
%
%   See also BIRADS_FEATS BOUND_FEATS ECHO_FEATS ORIENT_FEATS SHAPE_FEATS
%
%
%   References:
%   ----------
%   W. K. Moon, C. M. Lo, et al. "Quantitative ultrasound analysis for 
%   classification of BI-RADS category 3 breast masses," J Digit Imaging,
%   vol. 26, pp. 1091-1098, 2013.
%
%   W.-C. Shen, R.-F. Chang, W. K. Moon, Y.-H. Chou, C.-S. Huang, "Breast 
%   ultrasound computer-aided diagnosis using bi-rads features," Acad Radiol,
%   vol. 14, no. 8, pp. 928-939, 2007.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   MARGIN_FEATS Version 1.0 (Matlab R2014a Unix)
%   December 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = margin_feats(BW)
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
r = max(D(:))+1; % Radius
BW3 = xor(~BW2,C>=r);
BW3 = bwareaopen(BW3,2);
D2 = bwdist(~BW3);
BW3(D2<2) = 0;
% Count undulations without considering slighter undulations
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
x = [U A MUA];
feats = {'mMU','mMA','mMUA'};