% SOFTMAXNORM Sigmoidal normalization.
%   Y = SOFTMAXNORM(X) normalizes data in X (N samples-by-D features) with
%   the hyperbolic tangent function. The output data range is [-1,+1].
%   This kind of normalization is adequate to reduce the influence of extreme 
%   values or outliers in the data without removing them from the dataset.
%
%   [Y,MEAN,STD] = SOFTMAXNORM(X) returns the mean and standard deviation of
%   each column in X.
%
%   Y = SOFTMAXNORM(X,[MEAN;STD]) normalizes data in X with the given mean and 
%   standard deviation values.
%
%   Example:
%   -------
%   X1 = rand(10,2);
%   [Y1,m,s] = softmaxnorm(X1);
%   X2 = rand(10,2);
%   Y2 = softmaxnorm(X2,[m;s]);
%
%   See also MINMAXNORM STATNORM
%
%
%   Reference:
%   ---------
%   K. L. Priddy, P. E. Keller, Artificial Neural Networks: An Introduction.
%   Bellingham, WA: SPIE-The Int. Soc. Optical Eng., 2005.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   SOFTMAXNORM Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [G,m,s] = softmaxnorm(F,stats)
if nargin == 1
    m = mean(F,1);
    s = std(F,1)+1e-9;
elseif nargin == 2
    m = stats(1,:);
    s = stats(2,:);
end
FM = bsxfun(@minus,F,m);
DV = bsxfun(@rdivide,FM,s);
E  = exp(-DV);
G = (1-E)./(1+E); % hyperbolic tangent function. Range [-1,+1]
% G = 1./(1+E); % Logistic sigmoid function. Range [0,+1]