% BOOTSTRAP632AUC AUC estimation using the .632+ bootstrap method.
%   AUCh = BOOTSTRAP632AUC(X,Y) estimates the AUC using the .632+
%   bootstrap method from data in X (N samples-by-D features) and the 
%   targets in Y (N samples-by-1), where Yi=1 is a health sample and Yi=2 
%   is an unhealth sample (i=1,...,N). The LDA classifier is used to predict 
%   the out-of-bag data. The default number of bootstrap samples is set to 1000.
%   
%   AUCh = BOOTSTRAP632AUC(X,Y,B) estimates the AUC using B bootstrap samples.
%
%   Example:
%   -------
%   load('bcwd.mat');
%   X = softmaxnorm(X);         % Normalize data
%   AUCh = bootstrap632AUC(X,Y);% Estimate the AUC bootstrap   
%
%   See also BOOTSTRAP632ERR
%
%
% References:
% ----------
%   Berkman Sahiner, Heang-Ping Chan, and Lubomir Hadjiiski, "Classifier 
%   performance prediction for computer-aided diagnosis using a limited
%   dataset," Medical Physics, vol. 35, no.4, pp. 1559-1570, 2008.
%
%   W. Gomez, W. C. A. Pereira, A. F. C. Infantosi, "Improving classification
%   performance of breast lesions on ultrasonography," Pattern Recognition, 
%   vol. 48, pp. 1121-1132, 2015.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   BOOTSTRAP632AUC Version 1.0 (Matlab R2014a Unix)
%   June 2017
%   Copyright (c) 2017, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [AUC_hat,AUC_Xb_Xb0] = bootstrap632AUC(X,Y,B,opt)
if nargin < 4
    opt = 'soft';
end
if nargin < 3
    B = 1000;
    opt = 'soft';
end
Y = Y(:);       % force column vector
N = numel(Y);   % # observations
% Computing the resubstitution AUC value
Model = trainLDA(X,Y);
Out   = classifyLDA(X,Model);
if strcmpi(opt,'soft')
    AUC_X_X = AUC(Out.Scores,Y);
elseif strcmpi(opt,'hard')
    PP = classperf(Out.Labels,Y);
    AUC_X_X = PP.AUC;
end
% Buliding B boostrap samples
[iTr,iTt] = BootstrapSamples(N,B);
% Bootstrap samples AUC values
AUC_Xb_Xb0 = zeros(1,B);
for b = 1:B
    Model = trainLDA(X(iTr(:,b),:),Y(iTr(:,b),:));
    Out   = classifyLDA(X(iTt(:,b),:),Model);
    if strcmpi(opt,'soft')
        AUC_Xb_Xb0(b) = AUC(Out.Scores,Y(iTt(:,b),:));
    elseif strcmpi(opt,'hard')
        PP = classperf(Out.Labels,Y(iTt(:,b),:));
        AUC_Xb_Xb0(b) = PP.AUC;
    end
end
AUC_Xb_Xb0_p = max(0.5,AUC_Xb_Xb0); % Equation 7
% Overfitting rate, Equation 9
R_b = zeros(1,B);
R_b(AUC_Xb_Xb0<=0.5) = 1;
idx = (AUC_X_X > AUC_Xb_Xb0) & (AUC_Xb_Xb0 > 0.5);
R_b(idx) = (AUC_X_X - AUC_Xb_Xb0(idx))./(AUC_X_X - 0.5);
% Alpha weight, Equation 8
alpha_b = 0.632./(1-0.368*R_b);
% AUC estimated with 0.632+ bootstrap method, Equation 6
AUC_pre = (1-alpha_b)*AUC_X_X + alpha_b.*AUC_Xb_Xb0_p;
AUC_hat = mean(AUC_pre);
%***********************************************************
function [indBtr,indBtt] = BootstrapSamples(N,B)
indBtr = randi(N,N,B);
indBtt = false(N,B);
for b = 1:B
    indBtt(:,b) = oob(b,indBtr,N);
end
%***********************************************************
% out-of-bag data
function ibicont = oob(i,indBtr,N)
ids = ones(N,1);
Nib = accumarray(indBtr(:,i),ids,[N,1]); 
ibicont = Nib==0;