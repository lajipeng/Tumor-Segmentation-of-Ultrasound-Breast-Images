% EQUIVELLIPSE Morphological features based on equivalent ellipse.
%   [X,FEATS] = EQUIVELLIPSE(BW) computes eight features from the binary
%   blob of a breast lesion BW: orientation, elliptic-normalized circumference, 
%   elliptic-normalized skeleton, long axis to short axis ratio, proportional 
%   distance between edges, shape class, major axis length, and minor axis length.
%   X is a numeric vector with the feature values and FEATS is a cell vector with 
%   the name of the features in the same order as in X.
%
%   [X,FEATS] = EQUIVELLIPSE(BW,FEATURES) computes a specific set of 
%   morphological features according with FEATURES. The set of valid 
%   strings includes (case insensitive):
%
%       'all'   - Entire feature set.
%       'angle'	- Angle of the major axis of the equivalent ellipse.
%       'enc'   - Elliptic-normalized circumference.
%       'ens'   - Elliptic-normalized skeleton.
%       'ls'    - Long axis to short axis ratio.
%       'pdist' - Proportional distance between edges.
%       'shape' - Shape class feature.
%       'maxax' - Major axis length of the equivalent ellipse.
%       'minax' - Minor axis length of the equivalent ellipse.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = equivellipse(BW);
%   % Same as [x,feats] = equivellipse(BW,'all')
%
%   Example 2: Compute the angle
%   ----------------------------
%   load('BUS02.mat');   
%   [x,feats] = equivellipse(BW,'angle');
%
%   Example 3: Compute elliptic-normalized skeleton and shape class
%   ---------------------------------------------------------------
%   load('BUS02.mat');   
%   [x,feats] = equivellipse(BW,'ens','shape');
%
%   See also CONVHULLDIFF GEOMETRIC NRL
%
%
%   References:
%   ----------
%   C.-M. Chen, Y.-H. Chou, K.-C. Han, et al., "Breast lesions on sonograms: 
%   computer-aided diagnosis with nearly setting-independent features and 
%   artificial neural networks," Radiology, vol. 226, no. 2, pp. 504-514,2003.
%
%   W.-C. Shen, R.-F. Chang, W.K. Moon, Y.-H. Chou, C.-S. Huang, "Breast ultrasound
%   computer-aided diagnosis using bi-rads features," Academic Radioliology,
%   vol. 14, no. 8, pp. 928-939, 2007.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   EQUIVELLIPSE Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = equivellipse(BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'Angle';'ENC';'ENS';'LS';'PDist';'Shape';'MaxAx';'MinAx'};      
idxStats = RequestedInputs(officialStats,opts{:});
%*********************************************************************
% Elipse equivalente
% Calcula el area y el perimetro del objeto binario original
BW = double(BW);
Pbw = regionprops(BW,'Area','Centroid','Perimeter');
xc = Pbw.Centroid(1); yc = Pbw.Centroid(2);
A = Pbw.Area;
[y,x] = find(BW);
[xx,yy] = meshgrid(1:size(BW,2),1:size(BW,1));
% Calcula los momentos de segundo orden del objeto binario original
Sxx = (1/A)*sum((x-xc).^2);
Syy = (1/A)*sum((y-yc).^2);
Sxy = (1/A)*sum((x-xc).*(y-yc));
% Calcula los coeficientes de la ecuacion general de la elipse
coef = (1/(4*Sxx*Syy - Sxy^2))*[Syy -Sxy;-Sxy Sxx];
a = coef(1,1); b = coef(1,2); c = coef(2,2);
% Calcula la elipse eqivalente del objeto binario
E = (a*(xx-xc).^2 + 2*b*(xx-xc).*(yy-yc) + c*(yy-yc).^2) < 1;
% Calcula el contorno y el perimetro de la elipse equivalente
Pe = regionprops(E,'Perimeter','MinorAxisLength','MajorAxisLength');
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Orientacion
theta = 0.5*atan((2*Sxy)/(Sxx-Syy));
if (Sxx > Syy) && (theta > 0)
    theta = 1*theta;
elseif (Sxx > Syy) && (theta < 0)
    theta = -1*theta;
elseif (Sxx < Syy) && (theta > 0)
    theta = pi/2 - theta;
elseif (Sxx < Syy) && (theta < 0)
    theta = pi/2 - (-1*theta);
end
O = theta*180/pi;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ENC
ENC = Pbw.Perimeter/Pe.Perimeter;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ENS
BW_skl = bwmorph(BW,'skel',Inf);
ENS = sum(BW_skl(:))/Pe.Perimeter;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% LS
LS = Pe.MajorAxisLength/Pe.MinorAxisLength;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Proportinal distance
S1 = bwperim(BW);
S2 = bwperim(E);
avdis1 = averagedist(S1,S2);
avdis2 = averagedist(S2,S1);
PD = 100*((avdis1 + avdis2)/(2*sqrt(bwarea(E)/pi)));
% %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Shape class
% Parametrizacion de los contornos de BW y elipse
junk = bwboundaries(BW);
cBW  = junk{1};
yBW  = cBW(:,1); xBW = cBW(:,2);
junk = bwboundaries(E);
cE  = junk{1};
yE  = cE(:,1); xE = cE(:,2);
% Vectores unitarios de BW
rBW = [xBW-xc yBW-yc];
nBW = sqrt(sum(rBW.^2,2));
uBW = rBW./(repmat(nBW,1,2)+eps);
% Vectores unitarios de CH
rE = [xE-xc yE-yc];
nE = sqrt(sum(rE.^2,2));
uE = rE./(repmat(nE,1,2)+eps);
% Distancia entre vectores unitarios
D1 = dist(uBW,uE');
[~,ind] = min(D1,[],2);
% Correspondencia entre puntos de BW y puntos en CH con la orientacion
% mas proxima
mdE = cE(ind,:);
% Distancia Euclidiana
D2 = sqrt((cBW(:,1)-mdE(:,1)).^2+(cBW(:,2)-mdE(:,2)).^2);
SC = mean(D2);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
aux_x = [O ENC ENS LS PD SC Pe.MajorAxisLength Pe.MinorAxisLength];
aux_f = {'Orientation' 'ENC' 'ENS' 'L:S' 'PropDist' 'ShapeClass' 'MaxAxLe' 'MinAxLe'};
x = aux_x(idxStats);
feats = aux_f(idxStats);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function avdis = averagedist(cs,cr)
[lseg,cseg] = find(cs);
[lreal,creal] = find(cr);
[Lseg,Lreal] = meshgrid(lseg,lreal);
[Cseg,Creal] = meshgrid(cseg,creal);
dist = sqrt((Lseg-Lreal).^2+(Cseg-Creal).^2);
clear Lseg Lreal Cseg Creal
d = min(dist);
avdis = mean(d);