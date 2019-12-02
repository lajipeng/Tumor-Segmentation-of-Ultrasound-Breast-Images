% BIRADS_FEATS Compute BI-RADS features.
%   [X,FEAT] = BIRADS_FEATS(I,BW) computes BI-RADS features according to the
%   BI-RADS lexicon defined for breast ultrasound, where I is the gray-scale 
%   image containing a lesion and BW is the binary shape of the lesion:
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
%   Orientation             
%                           Angle of major axis of equivalent ellipse
%                           Depth-to-width ratio
%   
%   Margin
%                           Number of undulations (U)
%                           Angularity feature (A)
%                           Sum of U and A
%       
%   Boundary
%                           Abrupt interface at 10 pixels
%                           Abrupt interface at 25% of mass area
%                           Normalized radial gradient
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
%   Posterior behavior
%                           Minimum side difference
%                           Difference between mass and posterior intensities
%   
%   Example:
%   -------
%   load('BUS01.mat');   
%   [x,feat] = birads_feats(I,Smanual);
%
%   See also BOUND_FEATS ECHO_FEATS MARGIN_FEATS ORIENT_FEATS SHAPE_FEATS
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
%   BIRADS_FEATS Version 1.0 (Matlab R2014a Unix)
%   December 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feats] = birads_feats(I,BW)
[xshape,fshape] = shape_feats(BW);      % Shape
[xorie,forie] = orient_feats(BW);       % Orientation
[xmarg,fmarg] = margin_feats(BW);       % Margin
[xboun,fboun] = bound_feats(I,BW);      % Boundary
[xecho,fecho] = echo_feats(I,BW);       % Echo pattern
[xpab,fpab] = pab(I,BW);                % Posterior behavior
% Features
x = [xshape xorie xmarg xboun xecho xpab];
feats = [fshape forie fmarg fboun fecho fpab];