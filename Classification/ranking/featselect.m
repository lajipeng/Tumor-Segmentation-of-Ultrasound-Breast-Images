% FEATSELECT Feature selection based on MRMR criterion.
%   [MODEL,FEATS,STATS] = FEATSELECT(X,Y) selects a feature subset
%   with the best classification performance. The algorithm is performed
%   as follows:
%
%       1. Rank data in X (N samples-by-D features) with the MRMR
%       criterion. The target class Y (N samples-by-1) is '1' for health 
%       samples and '2' for unhealth samples. By default, the MRMR
%       criterion based on mutual information is used.
%       2. From ranked data, D subsets are evaluated, that is, the
%       first d features are taken until D features are considered. A LDA 
%       classifier is trained for each subset with d-features. In this
%       step a training set is considered.
%       3. The classification performance of a subset with d-features is
%       evaluated by considering a validation set.
%       4. The feature subset with the best classification performance is
%       used to train the final LDA model.
%
%   The outputs are:
%       1. MODEL: it is a structure with the trained LDA model.
%       2. FEATS: it is a structure with the ranked features (1-by-D vector) 
%       and the cut-off for the first best ranked features.
%       3. STATS: The mean and standard deviation of each column in X used
%       to normalize new data with softmax normalization.
%   
%   NOTE 1: The k-fold cross-validation method (with k=10) is used to 
%   validate a subset with d-features. The area under the ROC curve is used 
%   to measure the classification performance.
%
%   NOTE 2: If the Parallel Computing Toolbox is available, a parallel pool
%   is automatically open to speed-up the training procedure.
%   
%   Example:
%   -------
%   load('bcwd.mat');
%   ho = crossvalind('HoldOut',Y,0.2);  % Hold-out 80-20%
%   Xtr = X(ho,:);                      % Training data
%   Ytr = Y(ho,:);                      % Training labels
%   [Model,Feats,Stats] = featselect(Xtr,Ytr); % Feature selection model
%   Xtt = softmaxnorm(X(~ho,:),Stats);  % Test data normalization
%   Xtt = Xtt(:,Feats.Ranking);         % Rank test data
%   Xtt = Xtt(:,1:Feats.CutOff);        % Get feature subset
%   Ytt = Y(~ho,:);                     % Test labels
%   Out = classifyLDA(Xtt,Model);       % Evaluate LDA model
%   perf = classperf(Out.Labels,Ytt);   % Evaluate classification performance
%
%   See also MRMR_CORR MRMR_MI TRAINLDA CLASSIFYLDA
%
%
%   References:
%   ----------
%   H. Peng, F. Long, C. Ding, "Feature selection based on mutual information 
%   criteria of max-dependency, max-relevance, and min-redundancy," IEEE 
%   Trans. Pattern Anal. Mach. Intell., vol. 27, no. 8, pp. 1226-1238, 2005.
%
%   W. Gomez, W. C. A. Pereira, A. F. C. Infantosi, "Improving classification
%   performance of breast lesions on ultrasonography," Pattern Recognition, 
%   vol. 48, pp. 1121-1132, 2015.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   FEATSELECT Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [Model,Feats,Stats] = featselect(X,Y)
% Cross-validation and grid-search for SVM tuning
K  = 10; % 10-folds (modify it if necessary)
kf = crossvalind('KFold',Y,K); % Cross-validation partitions
% Normalization
D = size(X,2); % Dimensionality
[X,ave,sd] = softmaxnorm(X);
% Rank data with MRMR based on mutual information
idx = mrmr_mi(X,Y,2,'q');
X = X(:,idx);
% Check if Parallel Computing Toolbox
res = isToolboxAvailable('Parallel Computing Toolbox','warning');
% If it is installed uses parfor for speed-up training
perf = zeros(1,D);
% If it is installed uses parfor for speed-up training
if res
    n = feature('numCores');        % Number of physical cores
    s = isempty(gcp('nocreate'));   % Check if pool is already open
    if s
        parpool('local',n); % Open parallel pool
    end
    parfor d = 1:D
        perf(d) = eval_subset(X,Y,d,kf,K);
    end 
    if s
        delete(gcp); % Close parallel pool
    end
else
    for d = 1:D
        perf(d) = eval_subset(X,Y,d,kf,K);
    end
end
% Get best performance
[~,k] = max(perf);
% Train final model
Model = trainLDA(X(:,1:k),Y);
Feats.CutOff = k;
Feats.Ranking = idx;
Feats.Performance = perf;
Stats = [ave;sd];
%************************************************************************
function mperf = eval_subset(X,Y,d,kf,K)
kperf = zeros(1,K);
for k = 1:K
    ik = kf == k;   % Current fold
    Xtr = X(~ik,:); % Training set
    Ytr = Y(~ik,:);
    Xvd = X(ik,:);  % Validation set
    Yvd = Y(ik,:);
    Model = trainLDA(Xtr(:,1:d),Ytr); 
    Out = classifyLDA(Xvd(:,1:d),Model);
    kperf(k) = AUC(Out.Scores,Yvd);
end
mperf = mean(kperf);