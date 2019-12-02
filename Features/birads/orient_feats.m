% ORIENT_FEATS Compute BI-RADS orientation features.
%   [X,FEAT] = ORIENT_FEATS(BW) computes BI-RADS orientation features according 
%   to the BI-RADS lexicon defined for breast ultrasound, where BW is the binary 
%   shape of the lesion:
%   
%   BI-RADS feature         Quantitative feature
%   ---------------         ----------------------------
%   Orientation             
%                           Angle of major axis of equivalent ellipse
%                           Depth-to-width ratio
%   
%   Example:
%   -------
%   load('BUS01.mat');   
%   [x,feat] = orient_feats(Smanual);
%
%   See also BIRADS_FEATS BOUND_FEATS ECHO_FEATS MARGIN_FEATS SHAPE_FEATS
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
%   ORIENT_FEATS Version 1.0 (Matlab R2014a Unix)
%   December 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = orient_feats(BW)
BW = double(BW);
Pbw = regionprops(BW,'Area','Centroid','Perimeter');
xc = Pbw.Centroid(1); yc = Pbw.Centroid(2);
A = Pbw.Area;
[y,x] = find(BW);
% Calcula los momentos de segundo orden del objeto binario original
Sxx = (1/A)*sum((x-xc).^2);
Syy = (1/A)*sum((y-yc).^2);
Sxy = (1/A)*sum((x-xc).*(y-yc));
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
% Depth to width ratio
[yBW,xBW] = find(BW);
xBWmax = max(xBW); xBWmin = min(xBW);
yBWmax = max(yBW); yBWmin = min(yBW);
DWR = (yBWmax-yBWmin)/(xBWmax-xBWmin);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Features
x = [O DWR];
feats = {'oAngle','oDWR'};