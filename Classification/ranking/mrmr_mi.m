% MRMR_MI Minimal-redundancy-maximal-relevance (MRMR) criterion based on mutual information.
%   INDICES = MRMR_MI(X,Y,Q,OPTION) ranks the feature space by using the MRMR
%   criterion based on mutual information, where X is the feature matrix
%   (N samples-by-D features), Y are the class labels (N samples-by-1), Q
%   is the number of quantization levels to discretize continous variables,
%   and OPTION is a character, if 'q' performs quotient criterion and if
%   'd' performs difference criterion. INDICES is a vector of size 1-by-D
%   with the ranked indices of X.
%
%   NOTE: The functions "estpab" and "estmutualinfo" were compiled from the 
%   mutual information computation, conditional probability and entropy 
%   estimation for discrete/categorical variables package by Hanchuan Peng
%   at http://home.penglab.com/proj/mRMR/
%
%   Example:
%   -------
%   load('bcwd.mat');   
%   idx = mrmr_mi(X,Y,2,'q');
%   disp(idx);
%
%   See also FEATSELECT MRMR_CORR
%
%
%   Reference:
%   ---------
%   H. Peng, F. Long, C. Ding, "Feature selection based on mutual information 
%   criteria of max-dependency, max-relevance, and min-redundancy," IEEE 
%   Trans. Pattern Anal. Mach. Intell., vol. 27, no. 8, pp. 1226-1238, 2005.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   MRMR_MI Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function fea = mrmr_mi(x,f,Q,opt)
D = size(x,2);
if strcmpi(opt,'d')
    opt = 'minus';
elseif strcmpi(opt,'q')
    opt = 'rdivide';
else
    error('Only options d and q are accepted.');
end
[t,r] = RelRedMI(x,f,Q);
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
%*******************************************************************
function [F,R] = RelRedMI(X,Y,Q)
d = size(X,2);
QX = qDiscretize(X,Q);
F = zeros(d,1);
R = zeros(d,d);
for i = 1:d
    F(i) = mutualinfo(QX(:,i),Y);
    for j = 1:d
        R(i,j) = mutualinfo(QX(:,i),QX(:,j));
    end
end
%************************************************************************
function h = mutualinfo(vec1,vec2)
[p12, p1, p2] = estpab(vec1,vec2);
h = estmutualinfo(p12,p1,p2); 
%************************************************************************
function b = qDiscretize(a,d)
if nargin < 2
    d = 3;
end
[n,dim] = size(a);
b = zeros(n,dim);
for i = 1:dim
   b(:,i) = doDiscretize(a(:,i),d);
end
b = b+1;
%************************************************************************
function y_discretized= doDiscretize(y,d)
% discretize a vector
ys = sort(y);
y_discretized = y;
pos = ys(round(length(y)/d *(1:d)));
for j = 1:length(y)
    y_discretized(j) = sum(y(j)>pos);
end