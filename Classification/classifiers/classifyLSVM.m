% CLASSIFYLSVM Classify data with linear support vector machine.
%   OUT = CLASSIFLYSVM(X,MODEL) classifies the data in X (N samples-by-D 
%   features) by using a linear support vector machine (SVM), where
%   MODEL is a structure with the SVM parameters. OUT is a structure 
%   with both the scores and the class labels assigned to each sample in X.
%   If OUT.Labels=1 is a benign lesion,whereas if OUT.Labels=2 is a malignant 
%   lesion
%
%   NOTE 1: The function "predict" was compiled from the LIBSVM library
%   by Chih-Chung Chang and Chih-Jen Lin and the source codes are available
%   at https://github.com/cjlin1/liblinear/
%   
%   Example:
%   -------
%   load('bcwd.mat');
%   ho = crossvalind('HoldOut',Y,0.2);  % Hold-out 80-20%
%   [Xtr,m,s] = softmaxnorm(X(ho,:));   % Training data normalization
%   Xtt = softmaxnorm(X(~ho,:),[m;s]);  % Test data normalization
%   Ytr = Y(ho,:);                      % Training targets
%   Ytt = Y(~ho,:);                     % Test targets
%   Model = trainLSVM(Xtr,Ytr);          % Train SVM
%   Out = classifyLSVM(Xtt,Model);       % Test SVM
%   perf = classperf(Out.Labels,Ytt);   % Evaluate classification performance
%
%   See also AUC CLASSPERF CLASSIFYLDA CLASSIFYRBFN CLASSIFYSVM CLASSIFYRF
%
%
%   References:
%   ----------
%   Theodoridis S, Koutroumbas K. Pattern recognition. 4th edition. 
%   Burlington, MA: Academic Press 2009.
%
%   R.-E. Fan, P.-H. Chen, C.-J. Lin, "Working set selection using second 
%   order information for training SVM," Journal of Machine Learning Research,
%   vol. 6, pp. 1889-1918, 2005.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   CLASSIFYLSVM Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function Out = classifyLSVM(Xtt,Model)
[~,~,Sc] = predict(rand(size(Xtt,1),1), sparse(Xtt), Model,'-q'); 
Yp = Sc>=0;
Out.Scores = Sc;
Out.Labels = double(Yp)+1;