% DELINTUMOR Overlay tumor countour.
%   J = DELINTUMOR(I,BW) overlays the tumor countour of binary image BW 
%   over the gray-scale image I. The output image J is overlayed input
%   image.
%
%   J = DELINTUMOR(I,BW,OPT) displays the overlayed image I in a new figure 
%   if OPT=1, otherwise the image is not shown. By default OPT=1.
%
%   J = DELINTUMOR(I,BW,OPT,FNAME) saves the overlayed image in the file
%   FNAME.PNG, where FNAME is a string with the name of the output file.
%
%
%   Example 1:
%   ---------
%   load('BUS01.mat');
%   delintumor(I,Smanual);
%
%   Example 2:
%   ---------
%   load('BUS01.mat');
%   J = delintumor(I,Smanual,0,'overlay');
%
%   See also CROPROI

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   DELINTUMOR Version 2.0 (Matlab R2014a Unix)
%   June 2017
%   Copyright (c) 2017, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function I1 = delintumor(I0,BW,opt1,fname)

if nargin < 3
   opt1 = 1;
   fname = [];
end
if nargin < 4
   fname = [];
end

I0 = double(I0);
I0 = uint8(255*I0/max(I0(:)));
[M,N] = size(I0);
xy  = bwboundaries(BW,'noholes');
P   = false(M,N);
ind = sub2ind([M N],xy{1}(:,1),xy{1}(:,2));
P(ind) = 1;
S1 = imdilate(P,ones(7));
S2 = imdilate(P,ones(3));
I1 = I0;
I1(S1) = 255;
I1(S2) = 0;

if opt1
    figure;
    subplot 121; imshow(I0,'InitialMagnification', 'fit'); title('Original');
    subplot 122; imshow(I1,'InitialMagnification', 'fit'); title('Overlayed');
end

if ~isempty(fname)
    imwrite(I1,[fname '.png'],'png');
end