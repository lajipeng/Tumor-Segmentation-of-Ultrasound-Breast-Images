% MULTIRANKLET Multichannel decomposition by ranklet transform.
%   J = MULTIRANKLET(I,RESOLUTIONS,Q) performs the multichannel decomposition
%   of input image I by the ranklet transform, where RESOLUTIONS is a vector 
%   1xn of integers with the resolution values, and Q is the quantization 
%   level of output images. If Q is omitted, the output images are in the 
%   range [-1,+1]. J is a cell array where the columns represent the 
%   resolutions and the rows the orientations (1st: horizontal, 2nd: vertical and
%   3th: diagonal).
%
%   [J,M] = MULTIRANKLET(I,RESOLUTIONS,Q) all the orientations (Horizontal, 
%   Vertical, and Diagonal) at the same resolution are averaged; hence, if 
%   the size of the input image I is NxM then the output M is a hypermatrix 
%   of size MxNxRESOLUTIONS.
%
%   NOTE 1: For Windows, to run function "unifiedrankletPThreadWinMEX001.mewx64"
%   first copy the files "pthreadGC2.dll" and "pthreadVC2.dll" (provided in the
%   ZIP archive "...\C functions\Ranklet_source.zip") to the directory 
%   "C:\Windows\System32".
%
%   NOTE 2: The functions "unifiedrankletPThreadMEX001" (for MAC or LINUX)
%   and "unifiedrankletPThreadWinMEX001" (for Windows) were compiled 
%   from the source C code in the package "Ranklet_source.zip" located at the
%   "C functions" directory of this toolbox.
%
%   Example:
%   -------
%   I = imread('BUSreal.tif');
%   [J,M] = multiranklet(I,[2 4 8],32);
%   figure;
%   subplot 231; imshow(mat2gray(J{1,2})); title('Resolution 4/Horizontal');
%   subplot 232; imshow(mat2gray(J{2,2})); title('Resolution 4/Vertical');
%   subplot 233; imshow(mat2gray(J{3,2})); title('Resolution 4/Diagonal');
%   subplot 234; imshow(mat2gray(M(:,:,1))); title('Resolution 2');
%   subplot 235; imshow(mat2gray(M(:,:,2))); title('Resolution 4');
%   subplot 236; imshow(mat2gray(M(:,:,3))); title('Resolution 8');
%
%   See also MULTILOGGABOR PHASECONG QUANTIZATION
%
%
%   References:
%   ----------
%   M. Masotti, "A ranklet-based image representation for mass classification
%   in digital mammograms," Medical Physics, vol. 33, pp. 3951-3961, 2006.
%
%   M. Yang, W. Moon, Y. Wang, et al., "Robust texture analysis using 
%   multi-resolution gray-scale invariant features for breast sonographic 
%   tumor diagnosis," IEEE Trans. Med. Imaging, vol. 32, no. 12, 
%   pp. 2262-2273, 2013.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   MULTIRANKLET Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Arturo Rodriguez Cristerna, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [RT,T] = multiranklet(I,r,q)
str = computer; % Verify the kind of operative system for ranklet C code
I = uint8(I);
R = r(end);
Ri = numel(r);
T = zeros(size(I,1),size(I,2),Ri);
I = padarray(I,[R R],255,'both');
[Mi,Ni] = size(I);
RT = cell(3,Ri);
xc = R+1:Ni-R;
yc = R+1:Mi-R;
for i = 1:Ri
    % Ranklet transform
    if strcmpi(str,'MACI64')||strcmpi(str,'GLNXA64')
        [H,V,D] = unifiedrankletPThreadMEX001(I,r(i));
    elseif strcmpi(str,'PCWIN64')
        [H,V,D] = unifiedrankletPThreadWinMEX001(I,r(i));
    end
    M = (H+V+D)/3;
    if nargin == 3
        T(:,:,i) = quantization(M(yc,xc),q);
        RT{1,i}  = quantization(H(yc,xc),q);
        RT{2,i}  = quantization(V(yc,xc),q);
        RT{3,i}  = quantization(D(yc,xc),q);
    else
        T(:,:,i) = M(yc,xc);
        RT{1,i}  = H(yc,xc);
        RT{2,i}  = V(yc,xc);
        RT{3,i}  = D(yc,xc);
    end
end