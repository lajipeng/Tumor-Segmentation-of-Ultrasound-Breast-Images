% MRMR_CORR Minimal-redundancy-maximal-relevance (MRMR) criterion based on correlation.
%   INDICES = MRMR_CORR(X,Y,OPTION) ranks the feature space by using the MRMR
%   criterion based on correlation, where X is the feature matrix
%   (N samples-by-D features), Y are the class labels (N samples-by-1),
%   and OPTION is a character, if 'q' performs quotient criterion and if
%   'd' performs difference criterion. INDICES is a vector of size 1-by-D
%   with the ranked indices of X.
%
%   NOTE: The redundancy is measured by the Pearson correlation, whereas the
%   relevance is measured by the point biserial correlation. In the original
%   paper of Ding and Peng, the F-statistic is used to measure the relevance.
%
%   Example:
%   -------
%   load('bcwd.mat');   
%   idx = mrmr_corr(X,Y,'q');
%   disp(idx);
%
%   See also FEATSELECT MRMR_MI
%
%
%   Reference:
%   ---------
%   C. Ding, H. Peng, "Minimum redundancy feature selection from microarray 
%   gene expression data," Journal of Bioinformatics and Computational Biology,
%   vol. 3, no. 2, pp. 185-205, 2005.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   MRMR_CORR Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [fea,val] = mrmr_corr(x,f,opt)
%x = statnorm(x);
D = size(x,2);
if strcmpi(opt,'d')
    opt = 'minus';
elseif strcmpi(opt,'q')
    opt = 'rdivide';
else
    error('Only options d and q are accepted.');
end
[t,r] = RelRedCorr(x,f);
d = 1:D;
[val,fea] = max(t);
t(fea) = [];
d(fea) = [];
b = r(:,fea);
b(fea,:) = [];
while ~isempty(d)
    [aux,i] = max(feval(opt,t,mean(b,2)));
    fea = cat(2,fea,d(i));
    val = cat(2,val,aux);
    t(i) = [];
    d(i) = [];
    c = r(:,i);
    b(i,:) = [];
    b = cat(2,b,c(d));
end
%************************************************************************
function [F,R] = RelRedCorr(X,Y)
d = size(X,2);
C = unique(Y);
K = numel(C);
Fj = zeros(d,K);
for i = 1:d
    for j = 1:K
        Yj = Y==C(j);
        Fj(i,j) = pointbiserial(Yj,X(:,i));
    end
end
F = mean(abs(Fj),2); 
R = abs(corr(X,'type','Pearson'));
%*******************************************************************
function r = pointbiserial(d,x)
% Convert numeric values to logicals
d = logical(d);
% Length 
n = length(d);
% Lengths of groups 0 and 1
n1 = sum(d);
n0 = sum(~d);
% Mean of groups 0 and 1
x1 = mean(x(d));
x0 = mean(x(~d));
% SD of y
sx  = std(x);
%Correlation coefficient
r = (x1-x0)/sx*sqrt((n0*n1)/(n*(n-1)));