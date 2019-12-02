% CLASSIFYLDA Classify data with Fisher's linear discriminant analysis.
%   OUT = CLASSIFYLDA(X,MODEL) classifies the data in X (N samples-by-D 
%   features) by using a Fisher's linear discriminant analysis (LDA),
%   classifier where MODEL is a structure with the LDA parameters: weights 
%   and bias. OUT is a structure with both the scores and the class labels 
%   assigned to each sample in X. If OUT.Labels=1 is a benign lesion,
%   whereas if OUT.Labels=2 is a malignant lesion
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
%   perf = classperf(Out.Labels,Ytt);   % Evaluate classification performance
%
%   See also AUC CLASSPERF CLASSIFYRBFN CLASSIFYSVM CLASSIFYLSVM CLASSIFYRF
%
%
%   Reference:
%   ---------
%   Theodoridis S, Koutroumbas K. Pattern recognition. 4th edition. 
%   Burlington, MA: Academic Press 2009.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   CLASSIFYLDA Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function Out = classifyLDA(Xtt,Model)
W = Model.Weights;
b = Model.bias;
Sc = bsxfun(@minus,Xtt',b)'*W;
Yp = Sc>=0;
Out.Scores = Sc;
Out.Labels = double(Yp)+1;
