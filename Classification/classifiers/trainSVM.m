% TRAINSVM Train non-linear support vector machine classifier.
%   MODEL = TRAINSVM(X,Y) trains a support vector machine (SVM) with Gaussian
%   kernel from data in X (N samples-by-D features) and the targets  in Y 
%   (N samples-by-1), where Yi=1 is a health sample and Yi=2 is an unhealth
%   sample (i=1,...,N). MODEL a structure with the SVM model. The C and
%   gamma parameters are automatically adjusted by grid-search and 
%   cross-validation method.
%
%   [MODEL,PARAMS] = TRAINSVM(X,Y) obtains the C and gamma parameters found
%   by the grid-search procedure.
%
%   [MODEL,PARAMS,PERF] = TRAINSVM(X,Y) obtains the surface of AUC performance 
%   of all the points in the grid.
%
%   MODEL = TRAINSVM(X,Y,PARAMS) trains a support vector machine (SVM) with 
%   Gaussian kernel using the parameters C and gamma in the vector PARAMS,
%   such that, PARAMS = [C,gamma].
%
%   NOTE 1: The function "svmtrain" was compiled from the LIBSVM library
%   by Chih-Chung Chang and Chih-Jen Lin and the source codes are available
%   at http://www.csie.ntu.edu.tw/~cjlin/libsvm/
%   
%   NOTE 2: If PARAMS is not given, the SVM is trained with the k-fold 
%   cross-validation method (with k=10) and grid search to determine the  
%   parameters C and gamma. The search ranges are C = [2^-5,...,2^15] and 
%   gamma = [2^-15,...,2^3]. The area under the ROC curve is used to measure  
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
%   [Model,pars,S] = trainSVM(Xtr,Ytr); % Train SVM
%   Out = classifySVM(Xtt,Model);       % Test SVM
%   perf = classperf(Out.Labels,Ytt);   % Evaluate classification performance
%   figure; surf(-15:3,-5:15,S);        % Plot AUC surface
%   xlabel('log_2 \gamma');
%   ylabel('log_2 C');
%   zlabel('AUC value');
%   
%   Example 2:
%   ---------
%   load('bcwd.mat');
%   ho = crossvalind('HoldOut',Y,0.2);  % Hold-out 80-20%
%   [Xtr,m,s] = softmaxnorm(X(ho,:));   % Training data normalization
%   Xtt = softmaxnorm(X(~ho,:),[m;s]);  % Test data normalization
%   Ytr = Y(ho,:);                      % Training targets
%   Ytt = Y(~ho,:);                     % Test targets
%   Model = trainSVM(Xtr,Ytr,[64 0.125]); % Train SVM with C and gamma
%   Out = classifySVM(Xtt,Model);       % Test SVM
%   perf = classperf(Out.Labels,Ytt);   % Evaluate classification performance
%  
%   See also TRAINLDA TRAINRBFN TRAINLSVM TRAINRF
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
%   TRAINSVM Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [Model,params,mperf] = trainSVM(X,Y,params)
if nargin < 3
    % Cross-validation and grid-search for SVM tuning
    K  = 10; % 10-folds (modify it if necessary)
    kf = crossvalind('KFold',Y,K); % Cross-validation partitions
    % C and gamma ranges for grid-search
    % http://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf
    % (modify ranges if necessary)
    Cr = 2.^(-5:15);
    Gr = 2.^(-15:3);
    nC = numel(Cr);
    nG = numel(Gr);
    perf = zeros(nC,nG,K);
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
            perf(:,:,k) = grid_search(X,Y,kf,k,Cr,Gr);
        end 
        if s
            delete(gcp); % Close parallel pool
        end
    else
        for k = 1:K
            perf(:,:,k) = grid_search(X,Y,kf,k,Cr,Gr);
        end
    end
    % Get the configuration with the maximum classification performance
    mperf = mean(perf,3);   % Mean of k-folds results
    [~,k] = max(mperf(:));
    [k1,k2] = ind2sub([nC nG],k);
    C = Cr(k1); % C parameter
    g = Gr(k2); % Gamma parameter
    % Train final classication model
else
    C = params(1);
    g = params(2);
    mperf = [];
end
Model = SVM(X,Y,C,g);
params = [C,g];
%**********************************************************************
function perf = grid_search(X,Y,kf,k,Cr,Gr)
% Split data in training and test
ik = kf == k;   % Current fold
Xtr = X(~ik,:); % Training set
Ytr = Y(~ik,:);
Xvd = X(ik,:);  % Validation set
Yvd = Y(ik,:);
% Grid-search loops
nC = numel(Cr);
nG = numel(Gr);
perf = zeros(nC,nG);
for i = 1:nC
    for j = 1:nG 
        Model = SVM(Xtr,Ytr,Cr(i),Gr(j));   % Training
        Out = classifySVM(Xvd,Model);       % Validation
        perf(i,j) = AUC(Out.Scores,Yvd);    % AUC value
    end
end
%**********************************************************************
function Model = SVM(Xtr,Ytr,C,g)
Ytr(Ytr==1) = -1; % Health
Ytr(Ytr==2) = +1; % Unhealth
% Unbalanced classes correction
ni = sum(Ytr==-1);
nj = sum(Ytr==+1);
ci = (ni+nj)/(2*ni);
cj = (ni+nj)/(2*nj); 
% Parameters of SVM-lib
str = ['-q -s 0 -t 2 -c ' num2str(C) ' -g ' num2str(g) ...
       ' -w1 ' num2str(cj) ' -w-1 ' num2str(ci)];
Model = svmtrain(Ytr, sparse(Xtr), str);