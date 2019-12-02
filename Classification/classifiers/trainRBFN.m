% TRAINRBFN Train radial basis function network.
%   MODEL = TRAINRBFN(X,Y) trains a radial basis function network (RBFN)
%   from data in X (N samples-by-D features) and the targets in Y (N samples-by-1),
%   where Yi=1 is a health sample and Yi=2 is an unhealth sample (i=1,...,N). 
%   MODEL is a structure with the RBFN parameters: centroids, spreads, and weights.
%   The number of hidden units is automatically adjusted by grid-search and 
%   cross-validation method.
%
%   [MODEL,H] = TRAINRBFN(X,Y) obtains the number H of hidden units.
%
%   [MODEL,H,PERF] = TRAINRBFN(X,Y) obtains the curve of AUC performance 
%   of all the points in the grid.
%
%   MODEL = TRAINRBFN(X,Y,H) trains a radial basis function network (RBFN) with H 
%   hidden units.
%
%   NOTE 1: The RBFN is trained with the k-fold cross-validation method 
%   (with k=10) to determine the number of hidden units. The search range is
%   from 3 to 2*sqrt(N). The number of output units are limited to two
%   classes. The area under the ROC curve is used to measure the classification 
%   performance at each RBFN configuration.
%
%   NOTE 2: If the Parallel Computing Toolbox is available, a parallel pool
%   is automatically open to speed-up the training procedure.
%   
%   Example 1:
%   ---------
%   load('bcwd.mat');
%   ho = crossvalind('HoldOut',Y,0.2);  % Hold-out 80-20%
%   [Xtr,m,s] = softmaxnorm(X(ho,:));   % Training data normalization
%   Xtt = softmaxnorm(X(~ho,:),[m;s]);  % Test data normalization
%   Ytr = Y(ho,:);                      % Training targets
%   Ytt = Y(~ho,:);                     % Test targets
%   [Model,H,S] = trainRBFN(Xtr,Ytr);   % Train RBFN
%   Out = classifyRBFN(Xtt,Model);      % Test RBFN
%   perf = classperf(Out.Labels,Ytt);   % Evaluate classification performance
%   figure; plot(S,'k','linewidth',3);  % Plot AUC curve
%   xlabel('Number of hidden units');
%   ylabel('AUC value');
%
%   Example 2:
%   ---------
%   load('bcwd.mat');
%   ho = crossvalind('HoldOut',Y,0.2);  % Hold-out 80-20%
%   [Xtr,m,s] = softmaxnorm(X(ho,:));   % Training data normalization
%   Xtt = softmaxnorm(X(~ho,:),[m;s]);  % Test data normalization
%   Ytr = Y(ho,:);                      % Training targets
%   Ytt = Y(~ho,:);                     % Test targets
%   Model = trainRBFN(Xtr,Ytr,20);      % Train RBFN with 20 hidden units
%   Out = classifyRBFN(Xtt,Model);      % Test RBFN
%   perf = classperf(Out.Labels,Ytt);   % Evaluate classification performance
%
%   See also TRAINLDA TRAINSVM TRAINLSVM TRAINRF
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
%   TRAINRBFN Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [Model,H,mperf] = trainRBFN(X,Y,H)
if nargin < 3
    % Cross-validation and grid-search for RBFN hidden units tuning
    K  = 10; % 10-folds (modify it if necessary)
    kf = crossvalind('KFold',Y,K); % Cross-validation partitions
    % modify hidden units range if necessary
    hu = 3:2*fix(sqrt(numel(Y)));
    nh = numel(hu);
    perf = zeros(nh,K);
    % Check if Parallel Computing Toolbox
    res = isToolboxAvailable('Parallel Computing Toolbox','warning');
    % If it is installed uses parfor for speed-up training
    if res
        n = feature('numCores');        % Number of physical cores
        s = isempty(gcp('nocreate'));   % Check if pool is already open
        if s
            parpool('local',n); % Open parallel pool
        end
        parfor k = 1:K
            perf(:,k) = grid_search(X,Y,kf,k,hu);
        end
        if s
            delete(gcp); % Close parallel pool
        end
    else
        for k = 1:K
            perf(:,k) = grid_search(X,Y,kf,k,hu);
        end
    end
    % Get the configuration with the maximum classification performance
    mperf = mean(perf,2);   % Mean of k-folds results
    [~,k] = max(mperf);
    H = hu(k);
else
    mperf = [];
end
% Train final classication model
Model = RBFN(X,Y,H);
%**********************************************************************
function perf = grid_search(X,Y,kf,k,hu)
% Split data in training and test
ik = kf == k;   % Current fold
Xtr = X(~ik,:); % Training set
Ytr = Y(~ik,:);
Xvd = X(ik,:);  % Validation set
Yvd = Y(ik,:);
% Grid-search loop
nh = numel(hu);
perf = zeros(nh,1);
for i = 1:nh
    Model = RBFN(Xtr,Ytr,hu(i));   % Training
    if ~isempty(Model)
        Out = classifyRBFN(Xvd,Model); % Validation
        perf(i) = AUC(Out.Scores(:,2),Yvd); % AUC value to measure classification performance
    else
        perf(i) = 0;
    end
end
%*******************************************************************
function Model = RBFN(Xtr,Ytr,hh)
Ntr = numel(Ytr);
aux = double(Ytr==1);
bYtr = [aux 1-aux];
Si = zeros(1,hh);
t = 0;
while (~all(Si>0))&&(t<=20)
    t = t+1;
    Ci = mkMeans(Xtr,hh);
    DX = eucdist(Xtr,Ci);
    Si = get_sigmas(Ci);
end
if ~all(Si>0)
    Model = [];
    return;
end
S = Si(ones(Ntr,1),:);
aux1 = exp(-DX./(2*S.^2)); 
phi1 = cat(2,ones(Ntr,1),aux1);
W = pinv(phi1)*bYtr;
Model.Centroids = Ci;
Model.Sigmas = Si;
Model.Weights = W;
%*******************************************************************
% A Modified k-means Algorithm to Avoid Empty Clusters
% International Journal of Recent Trends in Engineering, Vol 1, No. 1, May 2009
function C = mkMeans(X,K)
n  = size(X,1);
ok = randperm(n,K);
C  = X(ok,:);
f  = 1;
tmax = 100;
t = 0;
while f && (t<=tmax)
    Cant = C;
    D = eucdist(X,C);
    [~,cl] = sort(D,2,'ascend');
    Y = cl(:,1);
    nk = accumarray(Y,ones(n,1),[K,1]);
    for i = 1:K
       C(i,:) = (1/(nk(i)+1))*(sum(X(Y==i,:),1) + C(i,:));
    end
    if norm(Cant-C) < 1e-9
        f = 0;
    end
    t = t+1;
end
%*******************************************************************
function S = get_sigmas(C)
DC = sqrt(eucdist(C,C));
DC = sort(DC,2,'ascend');
DC(:,1) = [];
S = mean(DC(:,1:2),2)';       % Mean distance to the 2-nearest neighbors
% S = sqrt(prod(DC(:,1:2),2))';   % Geometric mean distance to the 2-nearest neighbors
%*******************************************************************
function D = eucdist(A,B)
nA = sum(A.*A,2);
nB = sum(B.*B,2);
D = abs(bsxfun(@plus,nA,nB') - 2*(A*B'));