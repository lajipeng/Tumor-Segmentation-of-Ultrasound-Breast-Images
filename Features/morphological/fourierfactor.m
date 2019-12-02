% FOURIERFACTOR Fourier descriptors of lesion contour.
%   [X,FEAT] = FOURIERFACTOR(BW) computes the Fourier descriptors of the
%   binary blob of a breast lesion BW. X is the numeric value of the feature
%   and FEAT is the name of the feature.
%
%   Example:
%   -------
%   load('BUS02.mat');   
%   [x,feats] = fourierfactor(BW);
%
%   See also FOURIERSHAPE FRACTALCONTOUR POLYMODEL
%
%
%   Reference:
%   ---------
%   L. Shen, R. M. Rangayyan, J. L. Desautels, "Application of shape analysis
%   to mammographic calcifications," IEEE Trans. Med. Imaging, vol. 13, no. 2,
%   pp. 263-74, 1994.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   FOURIERFACTOR Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Arturo Rodriguez Cristerna, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [x,feat] = fourierfactor(BW)
Baux = bwboundaries(BW);
C = Baux{1};
z=C(:,2)+1i.*C(:,1);
N=size(C,1);
if mod(N,2)~=0
 z(N)=[];
 N=N-1;
end
Iu=zeros(1,N);

expValue=-2*1i*pi/N;
 for mu=1:N
     sumV=0;     
     for mx=1:N
         sumV=sumV+z(mx)*exp(expValue*(mu-1)*(mx-1));         
     end     
     Iu(1,mu)=sumV;
 end
 Iu=Iu./N; 
%Zo=abs(real(Iu));
Zo=abs((Iu));%absolute magnitudes 

Zo(1)=0;
Zo=Zo./Zo(2);
%trans
index1=1:floor(N/2);
index2=floor(N/2)+1:1:N;
Zto=[abs(Zo(index2)) abs(Zo(index1))];

divBase=[-size(index2,2)+1:1:-1 1 1:1:floor(N/2)];
zeroPos=size(index2,2);
if size(divBase,2)~=size(Zto,2)
    error('Error in the divbase calculus');
end
desc=abs(Zto)./(abs(divBase));
desc(zeroPos)=0;
sumDesc=sum(desc);
norm=sum(abs(Zto));
ff=(sumDesc./norm);%fourier factor
x=ff;
feat={'fourierFactor'};