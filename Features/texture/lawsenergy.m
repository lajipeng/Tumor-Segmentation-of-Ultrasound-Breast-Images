% LAWSENERGY Laws' texture energy measures.
%   [X,FEAT] = LAWSENERGY(I,BW,OPTION) computes five statistics of 15 
%   rotationally-invariant Laws' texture energy measures from the gray-level 
%   image I masked by the binary image BW. If OPTION=1 computes the features 
%   only from the internal region of the lesion and if OPTION=0 computes the 
%   features from the whole region of interest. X is a numeric vector with 
%   the feature values and FEATS is a cell vector with the name of the features 
%   in the same order as in X.
%
%   [X,FEATS] = LAWSENERGY(BW,FEATURES) computes a specific set of texture
%   features according with FEATURES. The set of valid strings includes 
%   (case insensitive):
%
%       'all'   - Entire feature set.
%       'mean'  - Mean value.
%       'std'   - Standard deviation.
%       'skw'   - Skew.
%       'kur'   - Kurtosis.
%       'egy'   - Energy.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS01.mat');   
%   [x,feats] = lawsenergy(I,Smanual,1);
%
%   Example 2: mean and standard deviation
%   --------------------------------------
%   load('BUS01.mat');   
%   [x,feats] = lawsenergy(I,Smanual,1,'mean','std');
%
%   See also AUTOCORR AUTOCOV BDIP_BVLC GLCM
%
%
%   References:
%   ---------
%   K. I. Laws, "Texture energy measures," In: DARPA Image Understanding
%   Workshop, Los Angeles, CA, 1979, pp. 47-51.
%
%   Kriti, J. Virmani, N. Dey, V. Kumar, "PCA-PNN and PCA-SVM Based CAD Systems 
%   for Breast Density Classification," in Applications of Intelligent Optimization
%   in Biology and Medicine: Current Trends and Open Problems, pp. 159-180,
%   2016.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   LAWSENERGY Version 2.0 (Matlab R2014a Unix)
%   December 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x2,feats2] = lawsenergy(I,BW,opt,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'mean','std','skw','kur','egy'};
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
opt = logical(opt);
I = double(I)+1;
[y,x] = find(BW);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
if exist('lawskernels.mat','file') == 2
    load('lawskernels.mat');
else
    K = lawskernels;
end

TIn = conv2(I,K{1,1},'same'); % Normalization mask L5L5 
% Compute texture energy images
TEM = cell(5,5);
h = ones(15);
for i = 1:5
   for j = 1:5 
       TI1 = conv2(I,K{i,j},'same')./TIn;
       TI2 = conv2(abs(TI1),h,'same');
       TEM{i,j} = TI2(ymin:ymax,xmin:xmax);
   end
end
% Compute rotationally-invariant texture energy images
[M,N] = size(TEM{1,1});
TR = zeros(M,N,15);
k = 0;
for i = 1:5
    for j = i:5
        k = k+1;
        TR(:,:,k) = 0.5*(TEM{i,j}+TEM{j,i});
    end
end
TR(:,:,1) = []; % Remove first image because it's always the unity
% Compute statistical parameters
fnames = {'eE5L5';'eS5L5';'eW5L5';'eR5L5';'eE5E5';'eS5E5';'eW5E5';'eR5E5';'eS5S5';'eW5S5';'eR5S5';'eW5W5';'eR5W5';'eR5R5'};
x = zeros(1,70);
feats = cell(1,70);
j = 1:5;
MN = M*N;
BW = BW(ymin:ymax,xmin:xmax);
for i = 1:14
    aux = TR(:,:,i);
    if opt
        TRi = aux(BW);
    else
        TRi = aux(:);
    end
    mn = sum(TRi)/MN;
    sd = sqrt(sum((TRi-mn).^2)/MN);
    x(j(1)) = mn;                                    feats(j(1)) = {[fnames{i} '_mean']};% Mean
    x(j(2)) = sd;                                    feats(j(2)) = {[fnames{i} '_std']}; % Standard deviation
    x(j(3)) = sum((TRi-mn).^3)/(MN*sd^3);            feats(j(3)) = {[fnames{i} '_skw']}; % Skewness
    x(j(4)) = (sum((TRi-mn).^4)/(MN*sd^4))-3;        feats(j(4)) = {[fnames{i} '_kur']}; % Kurtosis
    x(j(5)) = sum(TRi.^2)/MN;                        feats(j(5)) = {[fnames{i} '_egy']}; % Energy
    j = j+5;
end

feats2 = [];
x2 = [];
for i = 1:numel(idxStats)
    idx = not(cellfun(@isempty,strfind(feats,officialStats{idxStats(i)})'));
    feats2 = cat(2,feats2,feats(idx));
    x2 = cat(2,x2,x(idx));
end

%***************************************************************
function K = lawskernels

L5 = [1,4,6,4,1];   % Low-pass or average gray values
E5 = [-1,-2,0,2,1]; % Edges
S5 = [-1,0,2,0,-1]; % Spots
R5 = [1,-4,6,-4,1]; % Ripples
W5 = [-1,2,0,-2,1]; % Waves

K = cell(5,5);

K(1,1) = {L5'*L5};
K(2,1) = {L5'*E5};
K(3,1) = {L5'*S5};
K(4,1) = {L5'*W5};
K(5,1) = {L5'*R5};

K(1,2) = {E5'*L5};
K(2,2) = {E5'*E5};
K(3,2) = {E5'*S5};
K(4,2) = {E5'*W5};
K(5,2) = {E5'*R5};

K(1,3) = {S5'*L5};
K(2,3) = {S5'*E5};
K(3,3) = {S5'*S5};
K(4,3) = {S5'*W5};
K(5,3) = {S5'*R5};

K(1,4) = {W5'*L5};
K(2,4) = {W5'*E5};
K(3,4) = {W5'*S5};
K(4,4) = {W5'*W5};
K(5,4) = {W5'*R5};

K(1,5) = {R5'*L5};
K(2,5) = {R5'*E5};
K(3,5) = {R5'*S5};
K(4,5) = {R5'*W5};
K(5,5) = {R5'*R5};

save lawskernels.mat K;
