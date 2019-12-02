% NRL Morphological features based on the normalized radial lenght curve.
%   [X,FEATS] = NRL(BW) computes six morphological features based on the 
%   normalized radial lenght from the binary blob of a breast lesion BW: 
%   mean, standard deviation, area ratio, roughness, entropy, and zero-cross. 
%   X is a numeric vector with the feature values and FEATS is a cell vector 
%   with the name of the features in the same order as in X.
%
%   [X,FEATS] = NRL(BW,FEATURES) computes computes a specific set of
%   morphological features according with FEATURES. The set of valid 
%   strings includes (case insensitive): {'mean';'std';'ar';'rough','ent','zcr'};
%
%       'all'   - Entire feature set.
%       'mean'  - Mean value.
%       'std'   - Standard deviation.
%       'ar'    - Area ratio.
%       'rough' - Roughness.
%       'ent'   - Entropy.
%       'zcr'   - Zero cross.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = nrl(BW);
%   % Same as [x,feats] = nrl(BW,'all')
%
%   Example 2: Compute mean, standard deviation and entropy
%   -------------------------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = nrl(BW,'mean','std','ent');
%
%   See also CONVHULLDIFF EQUIVELLIPSE GEOMETRIC   
%
%
%   References:
%   ----------
%   A.V. Alvarenga, A.F.C. Infantosi, W.C.A. Pereira, C.M. Azevedo,
%   "Assessing the performance of morphological parameters in distinguishing 
%   breast tumors on ultrasound images," Med. Eng. Phys., vol. 32, no. 1, 
%   pp. 49-56, 2009.
%
%   Y.-H. Chou, C.-M. Tiu, G.-S. Hung, S.-C. Wu, T.Y. Chang, H.K. Chiang, 
%   "Stepwise logistic regression analysis of tumor contour features for 
%   breast ultrasound diagnosis," Ultrasound Med. Biol. vol. 27, pp. 1493-1498,
%   2001.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   MARGCLASS Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = nrl(BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'mean';'std';'ar';'rough';'ent';'zcr'};
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
BW = double(BW);
Prps = regionprops(BW,'Centroid');
xc = Prps.Centroid(1);
yc = Prps.Centroid(2);
Baux = bwboundaries(BW);
C = Baux{1};
% Normalized radial length
D = dist([xc yc],C');
Dn = D./max(D);
N = numel(Dn);
% Features
Dn_mean = mean(Dn); % Mean
Dn_sd   = std(Dn);  % Standard deviation
Dn_ar   = (1/(N*Dn_mean))*sum(Dn(Dn>Dn_mean)-Dn_mean); % Area ratio
Dn_rg   = (1/N)*sum(abs(diff(Dn))); % Roughness
h = histc(Dn,0:1/99:1);
Pk = h/sum(h);
Dn_ent  = -sum(Pk.*log(Pk+eps)); % Entropy 
Dn_zc = sum(abs(diff((Dn-Dn_mean)>0))); % Zero-cross
aux_x = [Dn_mean Dn_sd Dn_ar Dn_rg Dn_ent Dn_zc];
aux_f = {'sNRL_mean','sNRL_sd','sNRL_ar','sNRL_rg','sNRL_ent','sNRL_zc'};
x = aux_x(idxStats);
feats = aux_f(idxStats);