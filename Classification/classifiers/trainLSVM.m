% TRAINLSVM Train linear support vector machine classifier.
%   MODEL = TRAINLSVM(X,Y) trains a linear support vector machine (SVM)from 
%   data in X (N samples-by-D features) and the targets in Y (N samples-by-1),
%   where Yi=1 is a health sample and Yi=2 is an unhealth sample (i=1,...,N). 
%   MODEL a structure with the SVM model. The C parameter is automatically
%   adjusted by grid-search and cross-validation method.
%
%   [MODEL,C] = TRAINLSVM(X,Y) obtains the C parameter found by the
%   grid-search procedure.
%
%   [MODEL,PARAMS,PERF] = TRAINLSVM(X,Y) obtains the curve of AUC performance 
%   of all the points in the grid.
%
%   MODEL = TRAINSVM(X,Y,C) trains a linear support vector machine (SVM)
%   with the given C parameter.
%
%   NOTE 1: The function "train" was compiled from the LIBLINEAR library
%   by Chih-Chung Chang and Chih-Jen Lin and the source codes are available
%   at https://github.com/cjlin1/liblinear/
%   
%   NOTE 2: The SVM is trained with the k-fold cross-validation method 
%   (with k=10) and grid search to determine the parameter C in the range
%   C = [2^-5,...,2^15]. The area under the ROC curve is used to measure 
%   the classification performance at each SVM configuration.
%
%   NOTE 3: If the Parallel Computing Toolbox is available, a parallel pool
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
%   [Model,C,S] = trainLSVM(Xtr,Ytr);   % Train SVM
%   Out = classifyLSVM(Xtt,Model);      % Test SVM
%   perf = classperf(Out.Labels,Ytt);   % Evaluate classification performance
%   figure; plot(-5:15,S,'k','linewidth',3);    % Plot AUC curve
%   xlabel('log_2 C');
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
%   Model = trainLSVM(Xtr,Ytr,0.1250);  % Train SVM with C parameter
%   Out = classifyLSVM(Xtt,Model);      % Test SVM
%   perf = classperf(Out.Labels,Ytt);   % Evaluate classification performance
%
%   See also TRAINLDA TRAINRBFN TRAINSVM TRAINRF
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
%   TRAINLSVM Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [Model,C,mperf] = trainLSVM(X,Y,C)
if nargin < 3
    % Cross-validation and grid-search for SVM tuning
    K  = 10; % 10-folds (modify it if necessary)
    kf = crossvalind('KFold',Y,K); % Cross-validation partitions
    % C range for grid-search
    % http://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf
    % (modify ranges if necessary)
    Cr = 2.^(-5:15);
    nC = numel(Cr);
    perf = zeros(K,nC);
    % Check if Parallel Computing Toolbox
    res = isToolboxAvailable('Parallel Computing Toolbox','warning');
    % If it is installed uses parfor for speed-up training
    if res
        n = feature('numCores'); % Number of physical cores
        s = isempty(gcp('nocreate'));   % Check if pool is already open
        if s
            parpool('local',n); % Open parallel pool
        end
        parfor k = 1:K
            perf(k,:) = grid_search(X,Y,kf,k,Cr);
        end  
        if s
            delete(gcp); % Close parallel pool
        end
    else
        for k = 1:K
            perf(k,:) = grid_search(X,Y,kf,k,Cr);
        end
    end
    % Get the configuration with the maximum classification performance
    mperf = mean(perf,1);   % Mean of k-folds results
    [~,k] = max(mperf);
    C = Cr(k); % C parameter
else
    mperf = [];
end
% Train final classication model
Model = SVM(X,Y,C);
%**********************************************************************
function perf = grid_search(X,Y,kf,k,Cr)
% Split data in training and test
ik = kf == k;   % Current fold
Xtr = X(~ik,:); % Training set
Ytr = Y(~ik,:);
Xvd = X(ik,:);  % Validation set
Yvd = Y(ik,:);
% Grid-search loops
nC = numel(Cr);
perf = zeros(1,nC);
for i = 1:nC
    Model = SVM(Xtr,Ytr,Cr(i));   % Training
    Out = classifyLSVM(Xvd,Model);       % Validation
    perf(1,i) = AUC(Out.Scores,Yvd); % AUC value to measure classification performance
end
%**********************************************************************
function Model = SVM(Xtr,Ytr,C)
Ytr(Ytr==1) = -1; % Health
Ytr(Ytr==2) = +1; % Unhealth
% Unbalanced classes correction
ni = sum(Ytr==-1);
nj = sum(Ytr==+1);
ci = (ni+nj)/(2*ni);
cj = (ni+nj)/(2*nj); 
% Parameters of SVM-lib
str = ['-q -s 2 -c ' num2str(C) ' -w1 ' num2str(cj) ' -w-1 ' num2str(ci)];
Model = train(Ytr, sparse(Xtr), str);