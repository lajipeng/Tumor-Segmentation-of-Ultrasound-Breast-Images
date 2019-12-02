% BDIP_BVLC BDIP and BVLC moments
%   [X,FEATS] = BDIP_BVLC(I,BW) computes the BDIP and BVLC moments from
%   the gray-level image I masked by the binary image BW. X is a numeric 
%   vector with the feature values and FEATS is a cell vector with the name 
%   of the features in the same order as in X.
%
%   [X,FEATS] = BDIP_BVLC(BW,FEATURES) computes a specific set of texture
%   features according with FEATURES. The set of valid strings includes 
%   (case insensitive):
%
%       'all'   - Entire feature set.
%       'bdip'  - Difference of inverse probabilities.
%       'bvlc'  - Block variation of local correlation coefficients.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS01.mat');   
%   [x,feats] = bdip_bvlc(I,Smanual);
%   % Same as [x,feats] = bdip_bvlc(I,Smanual,'all')
%
%   Example 2: Compute BVLC moments
%   -------------------------------
%   load('BUS01.mat');   
%   [x,feats] = bdip_bvlc(I,Smanual,'bvlc');
%
%   See also AUTOCORR AUTOCOV GLCM LAWSENERGY  
%
%
%   References:
%   ----------
%   Y.-L. Huang, K.-L. Wang, D.-R. Chen, "Diagnosis of breast tumors with 
%   ultrasonic texture analysis using support vector machines," Neural Comput. 
%   Appl., vol. 15, no. 2, pp. 164-169, 2006.
%
%   Y. D. Chun, S. Y. Seo, N. C. Kim, "Image Retrieval Using BDIP and BVLC,"
%   IEEE Trans. Circuits Syst. Video Technol., vol. 13, no. 9, pp. 951-957,
%   2003.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   BDIP_BVLC Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = bdip_bvlc(US,BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'BDIP';'BVLC'};      
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
[y,x] = find(BW);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
I = double(US(ymin:ymax,xmin:xmax));
mn = min(I(:));
mx = max(I(:));
I = fix(255*((I-mn)/(mx-mn)));
M = 2; % Tamano del bloque
%------------------------------------------------------------------------- 
% Imagen BDIP
fun1 = @(block_struct) sum(sum(block_struct.data))*ones(size(block_struct.data));
fun2 = @(block_struct) max(max(block_struct.data))*ones(size(block_struct.data));
I1 = blockproc(I,[M M],fun1);
I2 = blockproc(I,[M M],fun2);
BDIP  = M.^2 - (I1./(I2+eps));
%------------------------------------------------------------------------- 
% Imagen BVLC
% Ventanas desplazadas
[Mi,Ni] = size(I);
yp = [1, 1:Mi-1];
ym = [2:Mi, Mi];
xp = [1, 1:Ni-1];
d10  = I(:,xp);
d01  = I(yp,:);
d11  = I(yp,xp);
d1_1 = I(ym,xp);
% Medias y desvios de las ventanas desplazadas
fun1  = @(block_struct) mean2(block_struct.data)*ones(size(block_struct.data));
fun2  = @(block_struct) std2(block_struct.data)*ones(size(block_struct.data));
mu10  = blockproc(d10,[M M],fun1);  sd10  = blockproc(d10,[M M],fun2);
mu01  = blockproc(d01,[M M],fun1);  sd01  = blockproc(d01,[M M],fun2);
mu11  = blockproc(d11,[M M],fun1);  sd11  = blockproc(d11,[M M],fun2);
mu1_1 = blockproc(d1_1,[M M],fun1); sd1_1 = blockproc(d1_1,[M M],fun2);
% Medias y desvios de la ventana original
mu00  = blockproc(I,[M M],fun1);  sd00  = blockproc(I,[M M],fun2);
% Computa VLCC
ps = zeros(Mi,Ni,4);
cte = 1/(M.^2);
fun1 = @(block_struct) sum(sum(block_struct.data))*ones(size(block_struct.data));
ps(:,:,1)  = (cte*blockproc((I.*d10 - mu00.*mu10),[M M],fun1))./(sd00.*sd10+eps);     % rho_(1,0)
ps(:,:,2)  = (cte*blockproc((I.*d01 - mu00.*mu01),[M M],fun1))./(sd00.*sd01+eps);     % rho_(0,1)
ps(:,:,3)  = (cte*blockproc((I.*d11 - mu00.*mu11),[M M],fun1))./(sd00.*sd11+eps);     % rho_(1,1)
ps(:,:,4)  = (cte*blockproc((I.*d1_1 - mu00.*mu1_1),[M M],fun1))./(sd00.*sd1_1+eps);  % rho_(1,-1)
% Computa BVLC
BVLC = max(ps,[],3) - min(ps,[],3);
%------------------------------------------------------------------------- 
% Clasificacion de los bloques basado en BDIP
[idC1,idC2,idC3,idC4,idC5,idC6,idC7,idC8] = classifying(BDIP);
% Momentos (m-media, s-desviacion estandar)
Dc = [mean(BDIP(idC1)) mean(BDIP(idC2)) mean(BDIP(idC3)) mean(BDIP(idC4)) mean(BDIP(idC5)) mean(BDIP(idC6)) mean(BDIP(idC7)) mean(BDIP(idC8)) ...
      std(BDIP(idC1))  std(BDIP(idC2))  std(BDIP(idC3))  std(BDIP(idC4))  std(BDIP(idC5))  std(BDIP(idC6))  std(BDIP(idC7))  std(BDIP(idC8))];
feaDc = {'BDIP_mC1' 'BDIP_mC2' 'BDIP_mC3' 'BDIP_mC4' 'BDIP_mC5' 'BDIP_mC6' 'BDIP_mC7' 'BDIP_mC8' ...
         'BDIP_sC1' 'BDIP_sC2' 'BDIP_sC3' 'BDIP_sC4' 'BDIP_sC5' 'BDIP_sC6' 'BDIP_sC7' 'BDIP_sC8'};  

Vc = [mean(BVLC(idC1)) mean(BVLC(idC2)) mean(BVLC(idC3)) mean(BVLC(idC4)) mean(BVLC(idC5)) mean(BVLC(idC6)) mean(BVLC(idC7)) mean(BVLC(idC8)) ...
      std(BVLC(idC1))  std(BVLC(idC2))  std(BVLC(idC3))  std(BVLC(idC4))  std(BVLC(idC5))  std(BVLC(idC6))  std(BVLC(idC7))  std(BVLC(idC8))];

feaVc = {'BVLC_mC1' 'BVLC_mC2' 'BVLC_mC3' 'BVLC_mC4' 'BVLC_mC5' 'BVLC_mC6' 'BVLC_mC7' 'BVLC_mC8' ...
         'BVLC_sC1' 'BVLC_sC2' 'BVLC_sC3' 'BVLC_sC4' 'BVLC_sC5' 'BVLC_sC6' 'BVLC_sC7' 'BVLC_sC8'};    

Dc(isnan(Dc))=0;%check for nan values
Vc(isnan(Vc))=0;

if numel(idxStats) == 2
    x = [Dc Vc];
    feats = [feaDc feaVc];
elseif idxStats == 1
    x = Dc;
    feats = feaDc;
elseif idxStats == 2
    x = Vc;
    feats = feaVc;    
end
%------------------------------------------------------------------------- 
%------------------------------------------------------------------------- 
function [idC1,idC2,idC3,idC4,idC5,idC6,idC7,idC8] = classifying(J)
% Clasificacion
% Primer paso
B1 = J >= mean2(J);
% Segundo paso
id1a = find(B1);
id2a = find(~B1);
mu1a = mean(J(id1a));
mu2a = mean(J(id2a));
id1b = id1a(J(id1a)>=mu1a);
id2b = id1a(J(id1a)<mu1a);
id3b = id2a(J(id2a)>=mu2a);
id4b = id2a(J(id2a)<mu2a);
% Tercer paso
muC12 = mean(J(id1b));
muC34 = mean(J(id2b));
muC56 = mean(J(id3b));
muC78 = mean(J(id4b));
idC1  = id1b(J(id1b)>=muC12); % Clase 1
idC2  = id1b(J(id1b)<muC12);
idC3  = id2b(J(id2b)>=muC34);
idC4  = id2b(J(id2b)<muC34);
idC5  = id3b(J(id3b)>=muC56);
idC6  = id3b(J(id3b)<muC56);
idC7  = id4b(J(id4b)>=muC78);
idC8  = id4b(J(id4b)<muC78);  % Clase 8