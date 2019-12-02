% STATNORM Statistical normalization.
%   Y = STATNORM(X) normalizes data in X (N samples-by-D features) to have 
%   zero mean and unit variance. This kind of normalization is also known 
%   as z-score normalization.
%
%   [Y,MEAN,STD] = STATNORM(X) returns the mean and standard deviation of
%   each column in X.
%
%   Y = STATNORM(X,[MEAN;STD]) normalizes data in X with the given mean and 
%   standard deviation values.
%
%   Example:
%   -------
%   X1 = rand(10,2);
%   [Y1,m,s] = statnorm(X1);
%   X2 = rand(10,2);
%   Y2 = statnorm(X2,[m;s]);
%
%   See also MINMAXNORM SOFTMAXNORM
%
%
%   Reference:
%   ---------
%   K. L. Priddy, P. E. Keller, Artificial Neural Networks: An Introduction.
%   Bellingham, WA: SPIE-The Int. Soc. Optical Eng., 2005.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   STATNORM Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [G,m,s] = statnorm(F,stats)
if nargin == 1
    m = mean(F,1);
    s = std(F,1)+1e-9;
elseif nargin == 2
    m = stats(1,:);
    s = stats(2,:);
end
FM = bsxfun(@minus,F,m);    % Zero mean
G  = bsxfun(@rdivide,FM,s); % Unit variance