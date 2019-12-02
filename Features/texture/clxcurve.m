% CLXCURVE Complexity curve.
%   [X,FEATS] = CLXCURVE(I,BW) computes five features from the complexity curve
%   of the gray-level image I masked by the binary image BW: Maximum value of 
%   transitions (mv), average value of transitions (av), aample mean (sm), 
%   sample standard deviation (ssd), and entropy (ent). X is a numeric vector 
%   with the feature values and FEATS is a cell vector with the name of the 
%   features in the same order as in X.
%
%   [X,FEATS] = CLXCURVE(BW,FEATURES) computes a specific set of texture
%   features according with FEATURES. The set of valid strings includes 
%   (case insensitive):
%
%       'all'   - Entire feature set.
%       'max'   - Maximum value of transitions.
%       'aver'  - Average value of transitions.
%       'smean' - Sample mean.
%       'std'   - Standard deviation.
%       'ent'   - Entropy.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS01.mat');   
%   [x,feats] = clxcurve(I,Smanual);
%   % Same as [x,feats] = clxcurve(I,Smanual,'all')
%
%   Example 2: Compute the entropy
%   ------------------------------
%   load('BUS01.mat');   
%   [x,feats] = clxcurve(I,Smanual,'ent');
%
%   Example 3: Compute the maximum, average, and standard deviation
%   ---------------------------------------------------------------
%   load('BUS01.mat');   
%   [x,feats] = clxcurve(I,Smanual,'max','aver','std');
%
%   See also HISTFEATURES FRACTALTEXTURE BGC
%
%   
%   Reference:
%   ---------
%   A.V. Alvarenga, W.C.A. Pereira, A.F.C. Infantosi, C.M. Azevedo, 
%   "Complexity curve and grey level co-occurrence matrix in the texture 
%   evaluation of breast tumor on ultrasound images," Medical Physics,
%   vol. 34, no. 2, pp. 379-387, 2007.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   CLXCURVE Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = clxcurve(US,BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'max';'aver';'smean';'std';'ent'}; 
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************

[y,x] = find(BW);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
I = double(US(ymin:ymax,xmin:xmax));
[Ny,Nx] = size(I);

ym = [2:Ny, Ny];
xm = [2:Nx, Nx];

grmn = min(I(:));
grmx = max(I(:));
I = fix(255*(I-grmn)/(grmx-grmn));
C = zeros(1,256);
cte_x = Nx*(Ny-1);
cte_y = Ny*(Nx-1);
cte   = cte_x + cte_y;
for i = 0:255
    B = double(I >= i);
    T10 = sum(sum(B-B(:,xm)<0));
    T01 = sum(sum(B-B(ym,:)<0));
    C(i+1) = (T10+T01)/cte;
end

mv = max(C);  % Maximum value of transitions
av = mean(C); % Average value of transitions
sm = (1/sum(C))*sum((0:255).*C); % Sample mean
ssd = sqrt((1/sum(C))*sum((((0:255)-sm).^2).*C)); % Sample standard deviation
ent = -sum(C(C>0).*log2(C(C>0))); % Entropy

aux_x = [mv av sm ssd ent];
aux_f = {'CC_mv' 'CC_av' 'CC_sm' 'CC_std' 'CC_ent'};
x = aux_x(idxStats);
feats = aux_f(idxStats);
