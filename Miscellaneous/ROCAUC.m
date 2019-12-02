% ROCAUC Receiver operating characteristic curve.
%   [Z,FPR,TPR] = ROCAUC(SCORES,LABELS) computes the receiver operating 
%   characteristic (ROC) curve and its corresponding area under the curve
%   (AUC) from classifier SCORES given the corresponding ground-truth LABELS. 
%   Z is the value of the area under the curve, FPR is the false posive rate
%   and TPR is the true positive rate.
%   
%   Example:
%   -------
%   load('bcwd.mat');
%   ho = crossvalind('HoldOut',Y,0.2);     % Hold-out 80-20%
%   [Xtr,m,s] = softmaxnorm(X(ho,:));      % Training data normalization
%   Xtt = softmaxnorm(X(~ho,:),[m;s]);     % Test data normalization
%   Ytr = Y(ho,:);                         % Training targets
%   Ytt = Y(~ho,:);                        % Test targets
%   Model = trainLDA(Xtr,Ytr);             % Train LDA
%   Out = classifyLDA(Xtt,Model);          % Test LDA
%   [z,FPR,TPR] = ROCAUC(Out.Scores,Ytt);  % ROC curve
%   figure;
%   plot(FPR,TPR,'k','linewidth',3); axis([0 1 0 1]);
%   title(['AUC = ' num2str(z)]);
%
%   See also CLASSPERF CLASSIFYLDA CLASSIFYRBFN CLASSIFYSVM CLASSIFYRF
%
%
%   Reference:
%   ---------
%   T. Fawcett, "An introduction to ROC analysis.," Pattern Recognition Letters,
%   vol. 27, no. 8, pp. 861-874, 2006. 

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   ROCAUC Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [croc,faRate,hitRate] = ROCAUC(score,trueLabel)
class1  = trueLabel==2;
Nclass1 = sum(class1);
class0  = trueLabel==1;
Nclass0 = sum(class0);
thresh  = unique(score);
thresh  = [thresh;max(score)+max(score)*0.01];
Nthresh = numel(thresh);
hitRate = zeros(1, Nthresh); 
faRate  = zeros(1, Nthresh);
for thi = 1:Nthresh
    th = thresh(thi);
    % hit rate = TP/P
    hitRate(thi) = sum(score(class1) >= th) / Nclass1;
    % fa rate = FP/N
    faRate(thi)  = sum(score(class0) >= th) / Nclass0;
end
croc = abs(trapz(faRate,hitRate));