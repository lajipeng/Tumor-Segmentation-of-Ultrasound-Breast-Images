% TRAINRF Train random forest classifier.
%   MODEL = TRAINRF(X,Y) trains a random forest classifier from data in 
%   X (N samples-by-D features) and the targets in Y (N samples-by-1), 
%   where Yi=1 is a health sample and Yi=2 is an unhealth sample (i=1,...,N). 
%   MODEL a structure with the random forest model. The default number of 
%   bootstrap samples is set to 1000.
%
%   MODEL = TRAINRF(X,Y,B) trains a random forest classifier with using B 
%   bootstrap samples.
%
%   NOTE 1: This function requires the "fitctree" function included in the 
%   Statistics and Machine Learning Toolbox.
%
%   NOTE 2: If the Parallel Computing Toolbox is available, a parallel pool
%   is automatically open to speed-up the training procedure.
%   
%   Example:
%   -------
%   load('bcwd.mat');
%   Model = trainRF(X,Y);            % Train random forest
%   Out = classifyRF(Model);         % Test random forest with out-of-bag data
%   perf = classperf(Out.Labels,Y);  % Evaluate classification performance
%  
%   See also TRAINLDA TRAINRBFN TRAINLSVM TRAINSVM
%
%
%   Reference:
%   ---------
%   Leo Breiman, "Random Forests", Machine Learning, vol. 45, no. 1, 
%   pp.  5-32, 2001.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   TRAINRF Version 1.0 (Matlab R2014a Unix)
%   June 2017
%   Copyright (c) 2017, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function Model = trainRF(X,Y,B)
if nargin < 3
    B = 1000;
end
[Data,Params] = split_bootstraps(X,Y,B);
B  = Params.B;
% Training with bootstrap samples
Forest = cell(1,B);
% Check if Parallel Computing Toolbox
res = isToolboxAvailable('Parallel Computing Toolbox','warning');
if res
    n = feature('numCores');        % Number of physical cores
    s = isempty(gcp('nocreate'));   % Check if pool is already open
    if s
        parpool('local',n); % Open parallel pool
    end
    parfor b = 1:B
        Forest{b} = train_trees(b,Data);
    end
    if s
        delete(gcp); % Close parallel pool
    end
else
    for b = 1:B
        Forest{b} = train_trees(b,Data);
    end    
end
Model.Forest = Forest;
Model.Data = Data;
Model.Params = Params;
% %***********************************************************************
function Models = train_trees(b,Data)
X  = Data.Train.X;
Y  = Data.Train.Y;
Ti = Data.Boost.Ti(:,b);
nfea = Data.nFeats;
Xtr = X(Ti,:);
Ytr = Y(Ti);
Models = train_tree(Xtr,Ytr,nfea);
%***********************************************************************
function Model = train_tree(X,Y,nfea)
Model = fitctree(X,Y,'MinLeaf',1,'Prune','off','Prior','uniform','NVarToSample',nfea);
%***********************************************************************
function [Data,Params] = split_bootstraps(Xtr,Ytr,B)
% Parameters
C = max(Ytr);       % Number of classes
[N,D] = size(Xtr);  % Number of patterns and features
% Parameters structure
Params.C = C;       % Number of classes
Params.D = D;       % Number of features
Params.N = N;       % Number of patterns
Params.B = B;       % Number of bootstraps
% Bootstrap indices
indBtr = randi(N,N,B);
ibicont = false(N,B);
for b = 1:B
    ibicont(:,b) = oob(b,indBtr,N);
end
% Number of random features
Feats = floor(sqrt(D)); 
% Feats = floor(1+log2(D)); % Original from Breiman's paper
% Data structures
Data.Train.X = Xtr;
Data.Train.Y = Ytr;
Data.Boost.Ti = indBtr;
Data.Boost.Vi = ibicont;
Data.nFeats = Feats;
%************************************************************************
% out-of-bag data
function ibicont = oob(i,indBtr,N)
ids = ones(N,1);
Nib = accumarray(indBtr(:,i),ids,[N,1]); 
ibicont = Nib==0;