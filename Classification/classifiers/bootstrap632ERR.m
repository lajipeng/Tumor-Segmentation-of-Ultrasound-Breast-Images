% BOOTSTRAP632ERR Classification error estimation using the .632+ bootstrap method.
%   E = BOOTSTRAP632ERR(X,Y) estimates the classification error using the 
%   .632+ bootstrap method from data in X (N samples-by-D features) and the 
%   targets in Y (N samples-by-1), where Yi=1 is a health sample and Yi=2 
%   is an unhealth sample (i=1,...,N). The LDA classifier is used to predict 
%   the out-of-bag data. The default number of bootstrap samples is set to 1000.
%   
%   E = BOOTSTRAP632ERR(X,Y,B) estimates the classification error using B 
%   bootstrap samples.
%
%   Example:
%   -------
%   load('bcwd.mat');
%   X = softmaxnorm(X);         % Normalize data
%   E = bootstrap632ERR(X,Y);	% Estimate the bootstrap error   
%
%   See also BOOTSTRAP632AUC
%
%
% References:
% ----------
%   B. Efron and R, Tibshirani, "Improvements on Cross-validation: 
%   the .632+ Bootstrap Method", Journal of the American Statistical 
%   Association, vol. 92, no. 438, pp. 548-560, 1997.
%
%   C. Ambroise and G. J. McLachlan, "Selection bias in gene extraction on
%   the basis of microarray gene-expression data", PNAS, vol. 99, no. 10,
%   pp.  6562-6566, 2002.
%
%   W. Gomez, W. Pereira, A.F.C. Infantosi, "Analysis of co-occurrence texture 
%   statistics as a function of gray-level quantization for classifying breast 
%   ultrasound," IEEE Trans. Med. Imaging, vol. 31, no. 10, pp. 1889-1899, 2012.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   BOOTSTRAP632ERR Version 1.0 (Matlab R2014a Unix)
%   June 2017
%   Copyright (c) 2017, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [e632p,err_boot] = bootstrap632ERR(X,Y,B)
if nargin < 3
    B = 1000;
end
Y = Y(:);       % force column vector
n = numel(Y);   % # observations
C  = unique(Y); % class labels
NC = numel(C);  % # classes
% Compute resubstitution error
Model = trainLDA(X,Y);
Out   = classifyLDA(X,Model);
Yrsb  = Out.Labels;
err  = Q(Y,Yrsb);
% Compute no-information rate
p     = zeros(1,NC);
q     = p;
gamma = 0;
for i = 1:NC % Ambroise & McLachlan, 2002
   p_idx = Y==C(i);      % input
   q_idx = Yrsb==C(i);   % output
   p(i)  = sum(p_idx)/n; % priors
   q(i)  = sum(q_idx)/n; % posteriors
   gamma = gamma + p(i)*(1-q(i)); % no-information rate
end
% Generate B bootstrap samples
[iTr,iTt] = genbootstrapsets(n,B);
% Compute leave-one-out bootstrap
err_boot = zeros(1,B);
for b = 1:B
    Model = trainLDA(X(iTr(:,b),:),Y(iTr(:,b),:));
    Out   = classifyLDA(X(iTt(:,b),:),Model);
    Ypred = Out.Labels;
    % Mis-classification rate
    err_boot(b) = Q(Y(iTt(:,b),:),Ypred);
end
% Bootstrap error
err_1 = mean(err_boot);
% Bootstrap .632
err_632 = (0.368*err) + (0.632*err_1);
% Relative overfitting rate
R = (err_1 - err)/(gamma - err);
% Make sure R ranges from 0-1
err_1 = min([err_1 gamma]);
if (err_1 > err) && (gamma > err)
   % R = R;
else
   R = 0;
end
% Bootstrap .632+
e632p = err_632 + (err_1 - err)*((0.368*0.632*R)/(1-0.368*R));
%************************************************************************
function [indBtr,indBtt] = genbootstrapsets(N,B)
indBtr = randi(N,N,B);
indBtt = false(N,B);
for b = 1:B
    indBtt(:,b) = oob(b,indBtr,N);
end
% out-of-bag
function ibicont = oob(i,indBtr,N)
ids = ones(N,1);
Nib = accumarray(indBtr(:,i),ids,[N,1]); 
ibicont = Nib==0;
%************************************************************************
% Discrepancy between actual and predicted labels
function r = Q(Ytrue,Ypred)
r = mean(Ytrue~=Ypred);