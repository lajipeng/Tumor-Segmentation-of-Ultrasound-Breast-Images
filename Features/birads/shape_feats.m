% SHAPE_FEATS Compute BI-RADS features.
%   [X,FEAT] = SHAPE_FEATS(BW) computes BI-RADS shape features according to the
%   BI-RADS lexicon defined for breast ultrasound, where BW is the binary shape 
%   of the lesion:
%   
%   BI-RADS feature         Quantitative feature
%   ---------------         ----------------------------
%   Shape                   
%                           Normalized residual value
%                           Normalized radial lenght
%                           Overlap with equivalent ellipse (EE)
%                           Elliptic-normalized circumference
%                           Elliptic-normalized skeleton
%                           Long axis to short axis ratio of EE
%                           Compactness or roundness
%                           Shape class
%                           Proportional distance between edges
%                           Major and minor axis length of EE
%
%   Example:
%   -------
%   load('BUS01.mat');   
%   [x,feat] = shape_feats(Smanual);
%
%   See also BIRADS_FEATS BOUND_FEATS ECHO_FEATS MARGIN_FEATS ORIENT_FEATS
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
%   SHAPE_FEATS Version 1.0 (Matlab R2014a Unix)
%   December 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = shape_feats(BW)
% Crop
[y,x] = find(BW);
xmn = min(x);
xmx = max(x);
ymn = min(y);
ymx = max(y);
BW2 = BW(ymn:ymx,xmn:xmx);
BW2 = padarray(BW2,[1 1],0,'both');
%---------------------------------------------------------------------
% Normalized Residual Value
BW_props = regionprops(BW2,'ConvexHull','Perimeter','Centroid','Area');
CH = roipoly(BW2,BW_props.ConvexHull(:,1),BW_props.ConvexHull(:,2));
NRV = bwarea(xor(BW2,CH))/bwarea(CH); 
%---------------------------------------------------------------------
% NRL
[xnrl,fnrl] = nrl(BW);
%---------------------------------------------------------------------
% Elipse equivalente
Pbw = regionprops(BW,'Area','Centroid','Perimeter');
A = Pbw.Area;
xc = Pbw.Centroid(1); yc = Pbw.Centroid(2);
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
E_props = regionprops(E,'Perimeter','MinorAxisLength','MajorAxisLength');
%---------------------------------------------------------------------
OR  = bwarea(BW&E)/bwarea(BW|E);	% Overlap
%---------------------------------------------------------------------
ENC = E_props.Perimeter/Pbw.Perimeter; % ENC
S = bwmorph(BW,'skel','inf');
S = bwareaopen(S,5);
ENS = sum(S(:))/E_props.Perimeter; % ENS
%---------------------------------------------------------------------
LS = E_props.MajorAxisLength/E_props.MinorAxisLength; %LS
Emax = E_props.MajorAxisLength;
Emin = E_props.MinorAxisLength;
%---------------------------------------------------------------------
R = 1 - ((4*pi*A)/(Pbw.Perimeter^2)); % Compactness
%---------------------------------------------------------------------
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
SC = sum(D2)/Pbw.Perimeter;
%---------------------------------------------------------------------
% Proportinal distance
S1 = bwperim(BW);
S2 = bwperim(E);
avdis1 = averagedist(S1,S2);
avdis2 = averagedist(S2,S1);
PD = 100*((avdis1 + avdis2)/(2*sqrt(bwarea(E)/pi)));
%---------------------------------------------------------------------
% Features
x = [NRV OR ENC ENS LS Emax Emin R SC PD xnrl];
feats = ['sNRV','sOR','sENC','sENS','sLS','sAX_MX','sAX_MN','sROUND','sSC','sPD',fnrl];
%---------------------------------------------------------------------
function avdis = averagedist(cs,cr)
[lseg,cseg] = find(cs);
[lreal,creal] = find(cr);
[Lseg,Lreal] = meshgrid(lseg,lreal);
[Cseg,Creal] = meshgrid(cseg,creal);
dist = sqrt((Lseg-Lreal).^2+(Cseg-Creal).^2);
clear Lseg Lreal Cseg Creal
d = min(dist);
avdis = mean(d);