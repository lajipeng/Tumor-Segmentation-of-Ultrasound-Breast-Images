% CLASSIFYRBFN Classify data with radial basis function network.
%   OUT = CLASSIFYRBFN(X,MODEL) classifies the data in X (N samples-by-D 
%   features) by using a radial basis function network (RBFN), where MODEL is 
%   a structure with the RBFN parameters: centroids, spreads, and weights. 
%   OUT is a structure with both the scores obtained from output units and the 
%   class labels assigned to each sample in X. If OUT.Labels=1 is a benign lesion,
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
%   Model = trainRBFN(Xtr,Ytr);         % Train RBFN
%   Out = classifyRBFN(Xtt,Model);      % Test RBFN
%   perf = classperf(Out.Labels,Ytt);   % Evaluate classification performance
%
%   See also AUC CLASSPERF CLASSIFYLDA CLASSIFYSVM CLASSIFYLSVM CLASSIFYRF
%
%
%   References:
%   ----------
%   Theodoridis S, Koutroumbas K. Pattern recognition. 4th edition. 
%   Burlington, MA: Academic Press 2009.
%
%   K. L. Priddy, P. E. Keller, Artificial Neural Networks: An Introduction.
%   Bellingham, WA: SPIE-The Int. Soc. Optical Eng., 2005.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   CLASSIFYRBFN Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function Out = classifyRBFN(Xtt,Model)
% Load model
Ci = Model.Centroids;
Si = Model.Sigmas;
W = Model.Weights;
% Classify
Ntt = size(Xtt,1);
DV = eucdist(Xtt,Ci);
S  = Si(ones(Ntt,1),:);
aux2 = exp(-DV./(2*S.^2));
phi2 = cat(2,ones(Ntt,1),aux2);
Out.Scores = phi2*W;
[~,Out.Labels] = max(Out.Scores,[],2);
%*******************************************************************
function D = eucdist(A,B)
nA = sum(A.*A,2);
nB = sum(B.*B,2);
D = abs(bsxfun(@plus,nA,nB') - 2*(A*B'));