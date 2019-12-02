% CLASSPERF Classification performance.
%   OUT = CLASSPERF(LPRED,LTRUE) measures distinct classification
%   performance indices, where LPRED are the predicted labels by the
%   classifier and LTRUE are the actual labels. OUT is a structure
%   containing the following indices: Matthews correlation coefficient (MCC),
%   accuracy (ACC), area under the ROC curve (AUC), sensitivity (SEN), 
%   specificity (SPE), positive predictive value (PPV), Negative predictive 
%   value (NPV), and confusion matrix.
%   
%   Example:
%   -------
%   load('bcwd.mat');
%   X = softmaxnorm(X);
%   Model = trainLDA(X,Y);           % Train LDA
%   Out = classifyLDA(X,Model);      % Test LDA
%   perf = classperf(Out.Labels,Y);  % Evaluate classification performance
%
%   See also AUC CLASSIFYLDA CLASSIFYRBFN CLASSIFYSVM CLASSIFYLSVM
%
%
%   Reference:
%   ---------
%   M. Sokolova, G. Lapalme, "A systematic analysis of performance measures
%   for classification tasks," Inform Process Manage, vol. 45, pp. 427-437,
%   2009.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   CLASSPERF Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function Out = classperf(Ypred,Ytrue)
pL = Ytrue == 2;
nL = Ytrue == 1;
TP  = sum(Ypred(pL)==2);
TN  = sum(Ypred(nL)==1);
FP  = sum(Ypred(nL)==2);
FN  = sum(Ypred(pL)==1);
num = (TP*TN)-(FP*FN);
den = sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
Out.MCC = max(0,num/(den+eps)); % Matthews correlation coefficient
Out.ACC = (TP+TN)/(TP+TN+FP+FN); % Accuracy                      
Out.AUC = 0.5*((TP/(TP+FN))+((TN/(TN+FP))));   % Area under the ROC curve
Out.SEN = TP/(TP+FN); % sensitivity or recall
Out.SPE = TN/(FP+TN); % specificity
Out.PPV = TP/(TP+FP); % precision or positive predictive value
Out.NPV = TN/(TN+FN); % negative prediction value
Out.Matrix = [TP FN;FP TN]; % Confusion matrix for binary classification