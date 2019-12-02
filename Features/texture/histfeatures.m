% HISTFEATURES First order statistics from gray-level histogram.
%   [X,FEATS] = HISTFEATURES(I,BW) computes 13 histogram features from the  
%   gray-level image I masked by the binary image BW: average histogram (AHg),
%   average gray level (AG), modified energy (MEgy), modified entropy (Metp), 
%   modified standard deviation (MSD), modified skew (MSk), average boundary 
%   gray level (BAG), difference (Diff), contrast (Ctr), energy (Egy), entropy
%   (Etp), standard deviation (SD), skew (Sk). X is a numeric vector with the feature
%   values and FEATS is a cell vector with the name of the features in the same 
%   order as in X.
%
%   [X,FEATS] = HISTFEATURES(BW,FEATURES) computes a specific set of texture
%   features according with FEATURES. The set of valid strings includes 
%   (case insensitive):
%   
%       'all'   - Entire feature set.
%       'bag'   - Average boundary grey level.
%       'ag'    - Average grey level.
%       'ahg'   - Average histogram.
%       'egy'   - Energy.
%       'megy'  - Modified energy.
%       'etp'   - Entropy.
%       'metp'  - Modified entropy.
%       'sd'    - Standard deviation.
%       'msd'   - Modified standard deviation.
%       'sk'    - Skew.
%       'msk'   - Modified skew.
%       'diff'  - Difference.
%       'ctr'   - Contrast.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS01.mat');   
%   [x,feats] = histfeatures(I,Smanual); 
%   % Same as [x,feats] = histfeatures(I,Smanual,'all')  
%
%   Example 2: Compute energy, entropy, and contrast
%   ------------------------------------------------
%   load('BUS01.mat');   
%   [x,feats] = histfeatures(I,Smanual,'egy','etp','ctr');   
%
%   See also CLXCURVE BGC FRACTALTEXTURE
%
%
%   Reference:
%   ---------
%   P. Zhang, B. Verma, K. Kumar, "Neural vs. statistical classifier in conjunction 
%   with genetic algorithm based feature selection," Pattern Recognit. Lett.,
%   vol. 26, no. 7, pp. 909-919, 2005.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   HISTFEATURES Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = histfeatures(US,BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'BAG';'AG';'AHg';'Egy';'MEgy';'Etp';'MEtp';'SD';'MSD';'Sk';'MSk';'Diff';'Ctr'};
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
US = double(US);
T = sum(BW(:)); % Numero total de pixeles en la lesion
K = 255;        % Numero total de numeles de gris
Ig = US(logical(BW)); % Niveles de gris en la lesion
C = imdilate(BW,ones(3))-imerode(BW,ones(3));
% Histograma solo de la region de la lesion
N = zeros(1,K+1);
for j = 0:K
    N(j+1) = sum(Ig==j);
end
Pj = N/T; % Probabilidad del nivel de gris
PIg = zeros(1,T);
for g = 1:T
    PIg(g) = N(Ig(g)+1)/T; 
end
j = 0:K;
BAG  = mean(US(logical(C)));
AG   = mean(Ig);
AHg  = mean(N);
Egy  = sum(Pj.^2);
MEgy = sum(PIg.^2);
Etp  = -sum(Pj(Pj>0).*log2(Pj(Pj>0)));
MEtp = -sum(PIg(PIg>0).*log2(PIg(PIg>0)));
SD   = sqrt(sum(((j-AG).^2).*Pj));
MSD  = sqrt(sum(((Ig-AG).^2).*PIg'));
Sk   = (1/(SD.^3))*sum((((j-AG).^3).*Pj));
MSk  = (1/(SD.^3))*sum(((Ig-AG).^3).*PIg');
Dffe = AG-BAG;
Ctr  = Dffe./(AG+BAG);
% Salidas
aux_x = [BAG AG AHg Egy MEgy Etp MEtp SD MSD Sk MSk Dffe Ctr];
aux_f = {'hBAG','hAG','hAHg','hEgy','hMEgy','hEtp','hMEtp','hSD','hMSD','hSk','hMSk','hDiff','hCtr'};
x = aux_x(idxStats);
feats = aux_f(idxStats);