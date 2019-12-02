% CLASSIFYRF Classify data with random forest classifier.
%   OUT = CLASSIFYRF(MODEL) classifies out-of-bag data by using a random
%   forest classifier, where MODEL is a structure with the RF parameters. 
%   OUT is a structure with both the votes and the class labels assigned to 
%   each sample in the original dataset. If OUT.Labels=1 is a benign lesion,
%   whereas if OUT.Labels=2 is a malignant lesion.
%
%   NOTE: This function requires the "predict" function included in the 
%   Statistics and Machine Learning Toolbox.
%   
%   Example:
%   -------
%   load('bcwd.mat');
%   Model = trainRF(X,Y);            % Train random forest
%   Out = classifyRF(Model);         % Test random forest with out-of-bag data
%   perf = classperf(Out.Labels,Y);  % Evaluate classification performance
%
%   See also AUC CLASSPERF CLASSIFYLDA CLASSIFYRBFN CLASSIFYLSVM CLASSIFYSVM
%
%
%   Reference:
%   ---------
%   Leo Breiman, "Random Forests", Machine Learning, vol. 45, no. 1, 
%   pp.  5-32, 2001.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   CLASSIFYRF Version 1.0 (Matlab R2014a Unix)
%   June 2017
%   Copyright (c) 2017, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function Out = classifyRF(Model)
B = Model.Params.B;
C = Model.Params.C;
X = Model.Data.Train.X;
Y = Model.Data.Train.Y;
N = numel(Y);
% classify out-of-bag
Ypb = zeros(N,B);
for i = 1:B
    Vi  = Model.Data.Boost.Vi(:,i);
    Xtt = X(Vi,:);
    Ypb(Vi,i) = predict(Model.Forest{i},Xtt);
end
% Majority vote
idx = not(logical(sum(abs(Ypb),2)));
Y(idx) = [];
Ypb(idx,:) = [];
N = numel(Y);
V = zeros(N,C);
for i = 1:N
    Ypi = Ypb(i,Model.Data.Boost.Vi(i,:));
    votes = accumarray(Ypi',ones(numel(Ypi),1),[C 1],@sum,0);
    V(i,:) = votes;
end
[~,Ypp] = max(V,[],2);
Out.Scores = V;
Out.Labels = Ypp;