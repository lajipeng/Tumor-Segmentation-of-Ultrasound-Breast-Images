% ECHO_FEATS Compute BI-RADS echo pattern features.
%   [X,FEAT] = ECHO_FEATS(I,BW) computes BI-RADS echo pattern features 
%   according to the BI-RADS lexicon defined for breast ultrasound, where 
%   I is the gray-scale image containing a lesion and BW is the binary shape 
%   of the lesion:
%   
%   BI-RADS feature         Quantitative feature
%   ---------------         ----------------------------
%   Echo pattern
%                           Average internal mass intensity
%                           Intensity difference between the 25 % 
%                               brighter pixels and whole tumor pixels
%                           GLCM features (D=[1,2,4]) from ranklets (r=[2,4,8]):
%                               Entropy, correlation, dissimilarity,
%                               homogeneity (mean, mean absolute deviation,
%                               and range over all orientations at same distance)
%                           Intensity autocorrelation
%                           Laws' energy measures features
%   
%   Example:
%   -------
%   load('BUS01.mat');   
%   [x,feat] = echo_feats(I,Smanual);
%
%   See also BIRADS_FEATS BOUND_FEATS MARGIN_FEATS ORIENT_FEATS SHAPE_FEATS
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
%   ECHO_FEATS Version 1.0 (Matlab R2014a Unix)
%   December 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = echo_feats(I,BW)
% Internal pattern
I = double(I);
J = I(BW)+1;
EPi = mean(J);
N = sum(BW(:));
H = accumarray(J,ones(N,1),[256 1],@sum,0);
p = 0.25*N;
for i = 256:-1:1
    if sum(H(i:256)) > p
        break;
    end
end
EPc = (mean(J(J>=i)) - EPi)/EPi;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% GLCM with ranklets
[~,J]  = multiranklet(I,[2 4 8]);
[y,x] = find(BW);
xmn = min(x);
xmx = max(x);
ymn = min(y);
ymx = max(y);
J2 = J(ymn:ymx,xmn:xmx,:);
BW2 = BW(ymn:ymx,xmn:xmx);
% Distance 1
[x11,f11] = glcm_feats(J2(:,:,1),BW2,1,'R2');
[x21,f21] = glcm_feats(J2(:,:,2),BW2,1,'R4');
[x31,f31] = glcm_feats(J2(:,:,3),BW2,1,'R8');
% Distance 2
[x12,f12] = glcm_feats(J2(:,:,1),BW2,2,'R2');
[x22,f22] = glcm_feats(J2(:,:,2),BW2,2,'R4');
[x32,f32] = glcm_feats(J2(:,:,3),BW2,2,'R8');
% Distance 4
[x14,f14] = glcm_feats(J2(:,:,1),BW2,4,'R2');
[x24,f24] = glcm_feats(J2(:,:,2),BW2,4,'R4');
[x34,f34] = glcm_feats(J2(:,:,3),BW2,4,'R8');
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Autocorrelacion
xa = autocorr(I,BW);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Laws energy
[xlaw,flaw] = lawsenergy(I,BW,1);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Features
x = [EPi EPc xa x11 x21 x31 x12 x22 x32 x14 x24 x34 xlaw];
feats = ['eEPi','eEPc','eACOR',f11,f21,f31,f12,f22,f32,f14,f24,f34,flaw];
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [x,feats] = glcm_feats(I,B,D,str)
warning off;
I(~B) = NaN;
Offset = D*[0 1; -1 1; -1 0; -1 -1]; % Distance 1
CMs = graycomatrix(I,'Of',Offset,'NumLevels',64,'GrayLimits',[-1 1]);
[size_CM_1,size_CM_2,~] = size(CMs);
[i,j] = meshgrid(1:size_CM_1,1:size_CM_2);
idx1 = (i+j)-1;
idx2 = abs(i-j)+1;
entro = zeros(1,4);
corrm = zeros(1,4);
dissi = zeros(1,4);
homom = zeros(1,4);
for l = 1:4 % number of CMs
    CMaux = CMs(:,:,l);
    % Normalize GLCM
    CM_sum = sum(CMaux(:));
    Pij = CMaux./CM_sum;  % Normalize each CM
    % 
    u_x = sum(sum(i.*Pij));
    u_y = sum(sum(j.*Pij));
    %
    p_xplusy  = zeros((2*size_CM_1 - 1),1); %[1]
    p_xminusy = zeros((size_CM_1),1);       %[1]
    for aux = 1:max(idx1(:))
       p_xplusy(aux) =  sum(Pij(idx1==aux));
    end
    for aux = 1:max(idx2(:))
       p_xminusy(aux) = sum(Pij(idx2==aux));
    end
    % Dissimilarity
    dissi(l) = sum(sum(abs(i-j).*Pij));
    % Homogeneity Matlab
    homom(l) = sum(sum(Pij./(1+abs(i-j)))); 
    % Entropy
    entro(l) = -sum(sum(Pij.*log(Pij+eps)));
    %
    s_x = sum(sum(Pij.*((i-u_x).^2)))^0.5;
    s_y = sum(sum(Pij.*((j-u_y).^2)))^0.5;
    % Correlation Matlab;
    corm = sum(sum(Pij.*(i-u_x).*(j-u_y)));
    corrm(l) = corm/((s_x*s_y)+0.0001);    
end
xm = [mean(entro) mean(corrm)...
      mean(dissi) mean(homom)]; % Average
xd = [mean(entro-xm(1)) mean(corrm-xm(2)) ...
      mean(dissi-xm(3)) mean(homom-xm(4))]; % Absolute deviation
xr = [max(entro)-min(entro) max(corrm)-min(corrm) ...
      max(dissi)-min(dissi) max(homom)-min(homom)]; % Range
str = ['D' num2str(D) '_' str];
featsm = {['eENTROm_' str],['eCORRm_' str],...
          ['eDISSm_' str],['eHOMOm_' str]};
featsd = {['eENTROd_' str],['eCORRd_' str],...
          ['eDISSd_' str],['eHOMOd_' str]};
featsr = {['eENTROr_' str],['eCORRr_' str],...
          ['eDISSr_' str],['eHOMOr_' str]};      
feats = [featsm featsd featsr];
x = [xm xd xr];