% TRAINLDA Train Fisher's linear discriminant analysis classifier.
%   MODEL = TRAINLDA(X,Y) trains a Fisher's linear discriminant analysis (LDA)
%   classifier from data in X (N samples-by-D features) and the targets in Y 
%   (N samples-by-1), where Yi=1 is a health sample and Yi=2 is an unhealth
%   sample (i=1,...,N). MODEL is a structure with the LDA parameters: weights 
%   and bias.
%   
%   Example:
%   -------
%   load('bcwd.mat');
%   ho = crossvalind('HoldOut',Y,0.2);  % Hold-out 80-20%
%   [Xtr,m,s] = softmaxnorm(X(ho,:));   % Training data normalization
%   Xtt = softmaxnorm(X(~ho,:),[m;s]);  % Test data normalization
%   Ytr = Y(ho,:);                      % Training targets
%   Ytt = Y(~ho,:);                     % Test targets
%   Model = trainLDA(Xtr,Ytr);          % Train LDA
%   Out = classifyLDA(Xtt,Model);       % Test LDA
%   perf = classperf(Out.Labels,Ytt);   % Evaluate classification performance
%
%   See also TRAINRBFN TRAINSVM TRAINLSVM TRAINRF
%
%
%   Reference:
%   ---------
%   Theodoridis S, Koutroumbas K. Pattern recognition. 4th edition. 
%   Burlington, MA: Academic Press 2009.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   TRAINLDA Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function Model = trainLDA(Xtr,Ytr)
Xtr = Xtr'; Ytr = Ytr';
n = size(Xtr,2);
mu1 = mean(Xtr(:,Ytr==1),2);    % Health
mu2 = mean(Xtr(:,Ytr==2),2);    % Unhealth
a1 = bsxfun(@minus,Xtr(:,Ytr==1),mu1);
a2 = bsxfun(@minus,Xtr(:,Ytr==2),mu2);
S_hat = (1/(n-2))*((a1*a1')+(a2*a2')); % Pooled covariance matrix
W = pinv(S_hat)*(mu2-mu1); % Weights
b = 0.5*(mu1+mu2); % Bias
Model.Weights = W;
Model.bias = b;