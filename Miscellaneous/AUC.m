% AUC Area under the ROC curve.
%   Z = AUC(SCORES,LABELS) computes the area under the ROC curve (AUC) from 
%   classifier SCORES given the corresponding ground-truth LABELS. Z is the 
%   value of the area under the curve.
%
%   Example:
%   -------
%   load('bcwd.mat');
%   ho = crossvalind('HoldOut',Y,0.2);  % Hold-out 80-20%
%   [Xtr,m,s] = softmaxnorm(X(ho,:));   % Training data normalization
%   Xtt = softmaxnorm(X(~ho,:),[m;s]);  % Test data normalization
%   Ytr = Y(ho,:);                      % Training targets
%   Ytt = Y(~ho,:);                     % Test targets
%   Model = trainLDA(Xtr,Ytr);          % Train LDA
%   Out = classifyLDA(Xtt,Model);       % Test LDA
%   Z = AUC(Out.Scores,Ytt);            % Measures de AUC
%
%   Reference:
%   ---------
%   D.J. Hand and R.J. Till, "A simple generalisation of the area under the 
%   ROC curve for multiple class classification problems," Machine Learning, 
%   vol. 45, no.2, pp. 171-186, 2001.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   AUC Version 1.0 (Matlab R2014a Unix)
%   June 2017
%   Copyright (c) 2017, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function out = AUC(Yp,Ytrue)
xi = Yp(Ytrue==2);
xj = Yp(Ytrue==1);
n0 = size(xi,1); 
n1 = size(xj,1);
r = tiedrank([xi;xj]);
s0 = sum(r(1:n0));
out = (s0-n0*(n0+1)/2)/(n0*n1);