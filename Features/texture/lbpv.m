% LBPV Local binary pattern variance from phase congruency.
%   [X,FEAT] = LBPV(I,BW) computes ten local binary pattern variance (lbpv) 
%   features from the phase congruency of the gray-level image I masked by  
%   the binary image BW. By default, six scales and eight orientations are 
%   used to compute the measure of phase congruency; hence, X is a vector 
%   of size 1x80, that is, ten lbpv features are computed for each phase 
%   congruency orientation. FEATS is a cell vector with the name of the features 
%   in the same order as in X.
%
%   [X,FEAT] = LBPV(I,BW,NSCALES,NORIENT) calculates the measure of phase 
%   congruency of the gray-scale image I with NSCALE scales and NORIENT
%   orientations; hence, X is a vector of size 1x(10xNORIENT).
%
%   Example:
%   -------
%   load('BUS01.mat');
%   [x,feats] = lbpv(I,Smanual);
%   
%   See also CLXCURVE BGC HISTFEATURES FRACTALTEXTURE
%
%
%   References:
%   ----------
%   L. Cai, X. Wang, Y. Wang, et al., "Robust phase-based texture descriptor 
%   for classification of breast ultrasound images," BioMedical Engineering 
%   OnLine, vol. 14, pp. 1-21, 2015.
%
%   Ojala T, Pietikainen M, Maenpaa T, "Multiresolution gray-scale and rotation 
%   invariant texture classification with local binary patterns," IEEE Trans 
%   Pattern Anal Mach Intell, vol. 24, no. 7, pp. 971-87, 2002.
%   
%   Guo Z, Zhang L, Zhang D, "Rotation invariant texture classification using 
%   LBP variance (LBPV) with global matching", Pattern Recognition, vol. 43
%   no. 3, pp. 706-719, 2010.
%
%   P. Kovesi, "Symmetry and asymmetry from local phase," In: Proceedings, 
%   10th Australian Joint Conference on Artificial Intelligence. 
%   Cambridge, MA: MIT Press; 1997, pp. 185-190.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   LBPV Version 1.0 (Matlab R2014a Unix)
%   June 2017
%   Copyright (c) 2017, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,f] = lbpv(I,Smanual,ns,no)
if nargin < 3
    ns = 6;
    no = 8;
end
% Crop ROI
J = cropROI(I,Smanual);
% Phase congruency
PC = phasecong(J,ns,no);
% Compute LBPV for each orientation
x = [];
f = [];
for i = 1:no
    [aux1,aux2] = PC_LBPV(PC{i},i);
    x = cat(2,x,aux1);
    f = cat(2,f,aux2);
end
%********************************************************************
function [xx,ff] = PC_LBPV(gc,o)
[Mi,Ni] = size(gc);
P = 8;
% P-neighbors
K = [0  0  0
     0  0  1
     0  0  0];
% Gray-Scale Invariance
s = zeros(Mi,Ni,P);
gp = zeros(Mi,Ni,P);
for p = 1:P
    gp(:,:,p) = imfilter(gc,K,'replicate');
    s(:,:,p) = double((gp(:,:,p)-gc)>=0);
    K = rot45(K);
end
% Uniform patterns
term1 = abs(s(:,:,P)-s(:,:,1));
term2 = zeros(Mi,Ni);
for p = 2:P
    term2 = term2 + abs(s(:,:,p)-s(:,:,p-1));
end
U = term1+term2;
% Rotation Invariance
aux = sum(s,3);
LBP_ = (P+1)*ones(Mi,Ni);
LBP_(U<=2) = aux(U<=2);
% Rotation Invariant Variance Measures
mu_ = mean(gp,3);
VAR_ = mean((gp-mu_(:,:,ones(1,P))).^2,3);
% LBP variance (LBPV)
[~,~,ks] = meshgrid(1:Ni,1:Mi,0:P+1);
idx = LBP_(:,:,ones(1,P+2))==ks;
xx = zeros(1,P+2);
ff = cell(1,P+2);
for k = 1:P+2
    aux = zeros(Mi,Ni);
    aux(idx(:,:,k)) = VAR_(idx(:,:,k));
    xx(k) = sum(aux(:));
    ff{k} = ['lbpv' num2str(k) '_o' num2str(o)];
end
%********************************************************************
function A = rot45(A)
idx = [2 1 4;3 5 7;6 9 8];
A = A(idx);