% BGC Binary gradient contours.
%   [X_BGC1,FEATS_BGC1,X_BGC2,FEATS_BGC2,X_BGC3,FEATS_BGC3,X_lbp,FEATS_LBP] = BGC(I,BW) 
%   computes four families of binary gradient contours (BGC) from the gray-level image 
%   I masked by the binary image BW: BGC1, BGC2, BGC3, and local binary pattern
%   (LBP). X is a numeric vector with the feature values and FEATS is a cell vector 
%   with the name of the features in the same order as in X.
%
%   [X,FEATS] = BGC(BW,FEATURES) computes a specific set of texture
%   features according with FEATURES. The set of valid strings includes 
%   (case insensitive):
%
%       'all'   - Entire feature set.
%       'bgc1'  - Binary gradient contours single-loop.
%       'bgc2'  - Binary gradient contours double-loop.
%       'bgc3'  - Binary gradient contours triple-loop.
%       'lbp'   - Local binary pattern.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS01.mat');   
%   [x_bgc1,feats_bgc1,x_bgc2,feats_bgc2,x_bgc3,feats_bgc3,x_lbp,feats_lbp] = bgc(I,Smanual);
%
%   Example 2: Compute bcg1
%   -----------------------
%   load('BUS01.mat');   
%   [x_bgc1,feats_bgc1] = bgc(I,Smanual,'bgc1');
%
%   Example 3: Compute bcg2 and lbp
%   -------------------------------
%   load('BUS01.mat');   
%   [x_bgc2,feats_bgc2,x_lbp,feats_lbp] = bgc(I,Smanual,'bgc2','lbp');
%
%   See also CLXCURVE HISTFEATURES LBPV FRACTALTEXTURE
%
%
%   Reference:
%   ---------
%   A. Fernandez, M. X. Alvarez, F. Bianconi, "Image classification with binary 
%   gradient contours," Optics and Lasers in Engineering, vol. 49,
%   pp. 1177-1184, 2011.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   BGC Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function varargout = bgc(I,BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'bgc1';'bgc2';'bgc3';'lbp'};
idxStats = RequestedInputs(officialStats,opts{:});
idxStats = sort(idxStats,'ascend');
%*********************************************************************
[y,x] = find(BW);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
I = double(I(ymin:ymax,xmin:xmax));
%I = double(I);
% Define shifted images
[Mi,Ni] = size(I);
S = [1,1:Mi-1];
N = [2:Mi, Mi];
E = [1,1:Ni-1];
W = [2:Ni, Ni];
J = zeros(Mi,Ni,8);
J(:,:,1) = I(:,W); % I0
J(:,:,2) = I(S,W); % I1
J(:,:,3) = I(S,:); % I2
J(:,:,4) = I(S,E); % I3
J(:,:,5) = I(:,E); % I4
J(:,:,6) = I(N,E); % I5
J(:,:,7) = I(N,:); % I6
J(:,:,8) = I(N,W); % I7

% Mask without image edges
idx = true(Mi,Ni);
idx(1,:)  = 0;
idx(:,1)  = 0;
idx(Mi,:) = 0;
idx(:,Ni) = 0;
Mi2 = Mi-2;
Ni2 = Ni-2;

xcell = cell(1,4);
fcell = cell(1,4);
if any(idxStats==1)
    % single-loop
    K = 2^8-1;
    BGC1_3x3 = zeros(Mi,Ni);
    for n = 0:7
        T1 = J(:,:,n+1);
        T2 = J(:,:,mod(n+1,8)+1);
        BGC1_3x3 = BGC1_3x3 + double(T1>=T2)*(2^n-1);
    end
    bgc1 = hist(BGC1_3x3(idx),0:K-1) / (Mi2*Ni2);
    feats_bgc1={};
    for m=0:254
        feats_bgc1={feats_bgc1{:} strcat('BGC1_',num2str(m))};
    end
    xcell{1} = bgc1;
    fcell{1} = feats_bgc1;
end

if any(idxStats==2)
    % double-loop
    K = (2^4-1)^2;
    BGC2_3x3a = zeros(Mi,Ni);
    BGC2_3x3b = zeros(Mi,Ni);
    for n = 0:3
        % patron rombo
        T1 = J(:,:,(2*n)+1);
        T2 = J(:,:,mod(2*(n+1),8)+1);
        BGC2_3x3a = BGC2_3x3a + double(T1>=T2)*2^n;
        % patron cuadrado
        T1 = J(:,:,(2*n+1)+1);
        T2 = J(:,:,mod(2*n+3,8)+1);
        BGC2_3x3b = BGC2_3x3b + double(T1>=T2)*(2^n-16);    
    end
    BGC2_3x3 = 15*BGC2_3x3a + BGC2_3x3b;
    bgc2 = hist(BGC2_3x3(idx),0:K-1) / (Mi2*Ni2);
    feats_bgc2={};
    for m=0:224
        feats_bgc2={feats_bgc2{:} strcat('BGC2_',num2str(m))};
    end
    xcell{2} = bgc2;
    fcell{2} = feats_bgc2;    
end

if any(idxStats==3)
    % triple-loop
    K = 2^8-1;
    BGC3_3x3 = zeros(Mi,Ni);
    for n = 0:7
        T1 = J(:,:,mod(3*n,8)+1);
        T2 = J(:,:,mod(3*(n+1),8)+1);
        BGC3_3x3 = BGC3_3x3 + double(T1>=T2)*(2^n-1);
    end
    bgc3 = hist(BGC3_3x3(idx),0:K-1) / (Mi2*Ni2);
    feats_bgc3={};
    for m=0:254
        feats_bgc3={feats_bgc3{:} strcat('BGC3_',num2str(m))};
    end
    xcell{3} = bgc3;
    fcell{3} = feats_bgc3;
end

if any(idxStats==4)
    % local binary pattern
    K = 2^8;
    LBP_3x3 = zeros(Mi,Ni);
    for n = 1:8
        LBP_3x3 = LBP_3x3 + double(J(:,:,n)>=I)*2^n;
    end
    lbp = hist(LBP_3x3(idx),0:K-1) / (Mi2*Ni2);
    feats_lbp={};
    for m=0:255
        feats_lbp={feats_lbp{:} strcat('LBP_',num2str(m))};
    end
    xcell{4} = lbp;
    fcell{4} = feats_lbp;
end


if numel(idxStats) == 4
    varargout{1} = bgc1;
    varargout{2} = feats_bgc1;
    varargout{3} = bgc2;
    varargout{4} = feats_bgc2;
    varargout{5} = bgc3;
    varargout{6} = feats_bgc3;     
    varargout{7} = lbp;
    varargout{8} = feats_lbp;
else
    non = 1;
    par = 2;
    for i = 1:numel(idxStats)
        varargout{non} = xcell{idxStats(i)};
        varargout{par} = fcell{idxStats(i)};
        non = non+2;
        par = par+2;
    end
end