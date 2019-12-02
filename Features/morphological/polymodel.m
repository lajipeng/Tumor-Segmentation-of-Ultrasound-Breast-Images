% POLYMODEL Polygonal model of lesion contour.
%   [X_TAF,FEATS_TAF,X_STAF,FEATS_STAF,X_PM,FEATS_PM] = POLYMODEL(BW) computes
%   three methods to analyze the  shape of a breast lesion BW: filtered 
%   turning angle function (TAF), smoothed turning angle function (STAF), and
%   polygonal model (PM). X is the numeric value of the features and FEATS is 
%   the name of the features in the same order as in X.
%
%   [X,FEATS] = POLYMODEL(BW,FEATURES) computes computes a specific set 
%   of morphological features according with FEATURES. The set of valid 
%   strings includes (case insensitive):
%
%       'all'   - Entire feature set.
%       'ftaf'  - Filtered turning angle function.
%       'staf'  - Smoothed turning angle function.
%       'poly'  - Polygonal model.
%
%   Example 1: Compute the entire feature set
%   -----------------------------------------
%   load('BUS02.mat');   
%   [x_TAF,feats_TAF,x_STAF,feats_STAF,x_PM,feats_PM] = polymodel(BW);
%   % Same as [x_TAF,feats_TAF,x_STAF,feats_STAF,x_PM,feats_PM] = polymodel(BW,'all')
%
%   Example 2: Compute the filtered turning angle function
%   ------------------------------------------------------
%   load('BUS02.mat');   
%   [x_TAF,feats_TAF] = polymodel(BW,'ftaf');
%
%   Example 3: Compute the smoothed turning angle function and polygonal model
%   --------------------------------------------------------------------------
%   load('BUS02.mat');   
%   [x_STAF,feats_STAF,x_PM,feats_PM] = polymodel(BW,'staf','poly');
%
%   See also FOURIERFACTOR FOURIERSHAPE FRACTALCONTOUR
%
%
%   Reference:
%   ---------
%   D. Guliato, R. M. Rangayyan, "Modeling and Analysis of Shape
%   with Applications in Computer-Aided Diagnosis of Breast Cancer,"
%   Morgan and Claypool Publishers, 2011.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   POLYMODEL Version 2.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Arturo Rodriguez Cristerna, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

%Polygonal Modeling of Contours of Breast Tumors With the Preservation of Spicules
%Denise Guliato, Rangaraj M. Rangayyan, Fellow, IEEE, Juliano D. Carvalho, and S??rgio A. Santiago
%It was added the smoothed TAF (STAF) descriptors of 
%Modeling shape and analysis of shape breast cancer - Guliato and rangaraj - 2011.pdf
%fTAF - filtered TAF
%STAF - smoothed TAF
%PM - polygonal model
function varargout = polymodel(BW,varargin)
%*********************************************************************
% Check for options
if nargin < 2
    opts = {'all'};
else
    opts = varargin;
end
officialStats = {'ftaf';'staf';'poly'};
idxStats = RequestedInputs(officialStats,opts{:});
idxStats = sort(idxStats,'ascend');
%*********************************************************************
% Baux = bwboundaries(BW);
% C = Baux{1};
[BWy, BWx]=size(BW);
for row = BWy:-1:1
    for col = BWx:-1:1
        if BW(row,col)==1
            break;
        end
    end
    if BW(row,col)==1
        break;
    end
end
C = bwtraceboundary(BW,[row col],'W',8,inf,'counterclockwise');
[BWy, BWx]=size(BW);
r = max_inner_circle(C(:,2),C(:,1));
%%display('Radio of the maximum inner circle,',num2str(r));
np=size(C,1);
id1=[1:np];
id2=[1,1:np-1];
%turning pixel
d(:,1)=C(id1,1)-C(id2,1);%y
d(:,2)=C(id1,2)-C(id2,2);%x
%turning angle
dtheta=rad2deg(atan2(d(:,1),d(:,2)));
ddtheta=dtheta(id1)-dtheta(id2);

pf.y=C(1,1); pf.x=C(1,2);
npplg=0;%Getting the line segments of the original contour
for m=1:np %each time that the direction of dtheta changes we add a new line segment
    if ddtheta(m)~=0
        line.pf=pf;%first point of the segment
        pl.y=C(m,1); pl.x=C(m,2);
        line.pl=pl;%last point of the segment
        line.dis= sqrt( (line.pl.x-line.pf.x)^2 + (line.pl.y-line.pf.y)^2 );%distance between points
        line.v1=line.pl.x-line.pf.x; %direction vector
        line.v2=line.pl.y-line.pf.y;
        npplg=npplg+1;
        plg{npplg}=line;
        pf=pl;%set the new first point
    end
end
%if not exists, then 
%we must add a line to the first point
if pf.y~=C(1,1)|| pf.x~=C(1,2)    
    line.pf=pf;%first point of the segment
    pl.y=C(1,1); pl.x=C(1,2);%returning to the origin
    line.pl=pl;%last point of the segment
    line.dis= sqrt( (line.pl.x-line.pf.x)^2 + (line.pl.y-line.pf.y)^2 );%distance between points
    line.v1=line.pl.x-line.pf.x; %direction vector
    line.v2=line.pl.y-line.pf.y;
    npplg=npplg+1;
    plg{npplg}=line;
end

%Set of distances and angles to use to obtain the filtered TAF
setOfDistances=[0,1,3,6]; %minDis=r/6;
setOfTheta=[150,160,170]; %maxTheta=170;
xTAFn=0; %number of times that we add descriptors
featsfTAF={}; xfTAF=[]; %features from the filtered TAF
featsSTAF={}; xSTAF=[]; %features of the STAF
featsPM={}; xPM=[]; %features of the polygonal model
%Get the filtered TAF
base_plg=plg;%storing the polygon obtained from the original contour
%Polygon of the original contour
Co=[];
Co(1,1)=base_plg{1}.pf.y; Co(1,2)=base_plg{1}.pf.x;%first point of the segment
for n=1:size(base_plg,2) %Get the points of the polygon
    Co(n+1,1)=base_plg{n}.pl.y; Co(n+1,2)=base_plg{n}.pl.x;%last point of the segment
end
BWCo = poly2mask(Co(:,2), Co(:,1), BWy, BWx);
for sod=1:size(setOfDistances,2);
    for sot=1:size(setOfTheta,2)
        plg={};
        plg=base_plg; %starting the computation of the TAF with the polygon obtained from the original contour
        if setOfDistances(sod)>0
           minDis=r/setOfDistances(sod); %The minumum distance could be based in the radio of the maxium inner circle - this was found manually
        else
           minDis=10;
        end
        maxTheta=setOfTheta(sot);
        olimit=size(plg,2)-1; limit=size(plg,2)-2;%just to start the loop
        while (olimit>limit)
            m=1; olimit=size(plg,2)-1; limit=size(plg,2)-1;
            while(m<limit)
                %First rule, if the two segments are shorter than a minimum
                %distance, then we join both segments
                if (plg{m}.dis<minDis) && (plg{m+1}.dis<minDis)
                    %if (plg{m}.dis+plg{m+1}.dis)<minDis
                    plg{m}.pl=plg{m+1}.pl;
                    plg{m}.dis= sqrt( (plg{m}.pl.x-plg{m}.pf.x)^2 + (plg{m}.pl.y-plg{m}.pf.y)^2 );
                    plg{m}.v1=plg{m}.pl.x-plg{m}.pf.x;
                    plg{m}.v2=plg{m}.pl.y-plg{m}.pf.y;
                    plg(m+1)=[];
                    limit=limit-1;
                else
                    %second rule, if one of the segments are longer than a minimum
                    %distance, then we analize the angle between them
                    %if the minim inner angle is higger than a maximum angle,
                    %then we join both segments
                    %angle between adjacent segments
                    %%%cosTheta=abs(plg{m}.v1*plg{m+1}.v1+plg{m}.v2*plg{m+1}.v2)/(sqrt(plg{m}.v1^2+plg{m}.v2^2)+sqrt(plg{m+1}.v1^2+plg{m+1}.v2^2)) esta forma no sirve error;
                    %http://www.jorge-fernandez.es/proyectos/angulo/temas/temap/index.html
                    a= sqrt((plg{m}.pf.x-plg{m}.pl.x)^2+(plg{m}.pf.y-plg{m}.pl.y)^2);
                    b= sqrt((plg{m+1}.pf.x-plg{m+1}.pl.x)^2+(plg{m+1}.pf.y-plg{m+1}.pl.y)^2);
                    c= sqrt((plg{m}.pf.x-plg{m+1}.pl.x)^2+(plg{m}.pf.y-plg{m+1}.pl.y)^2);
                    cosTheta= (c^2-b^2-a^2)./(-2*a*b);
                    theta=rad2deg(acos(cosTheta));
                    if (plg{m}.dis>minDis) || (plg{m+1}.dis>minDis)
                        if theta>maxTheta
                            plg{m}.pl=plg{m+1}.pl;
                            plg{m}.dis= sqrt( (plg{m}.pl.x-plg{m}.pf.x)^2 + (plg{m}.pl.y-plg{m}.pf.y)^2 );
                            plg{m}.v1=plg{m}.pl.x-plg{m}.pf.x;
                            plg{m}.v2=plg{m}.pl.y-plg{m}.pf.y;
                            plg(m+1)=[];
                            limit=limit-1;
                        end
                    end
                end
                m=m+1;%move forward in the loop
            end
            %this step is added to compute the "joining rules" between the
            %last and the firts segment
            if (plg{1}.dis<minDis) && (plg{end}.dis<minDis)
                %if (plg{1}.dis+plg{end}.dis)<minDis
                plg{1}.pf=plg{end}.pf;%change the first point
                plg{1}.dis= sqrt( (plg{1}.pl.x-plg{1}.pf.x)^2 + (plg{1}.pl.y-plg{1}.pf.y)^2 );
                plg{1}.v1=plg{1}.pl.x-plg{1}.pf.x;
                plg{1}.v2=plg{1}.pl.y-plg{1}.pf.y;
                plg(end)=[];
                limit=limit-1;
            else
                %%%cosTheta=abs(plg{1}.v1*plg{end}.v1+plg{1}.v2*plg{end}.v2)/(sqrt(plg{1}.v1^2+plg{1}.v2^2)+sqrt(plg{end}.v1^2+plg{end}.v2^2)) esta forma no sirve error;
                a= sqrt((plg{end}.pf.x-plg{end}.pl.x)^2+(plg{end}.pf.y-plg{end}.pl.y)^2);
                b= sqrt((plg{1}.pf.x-plg{1}.pl.x)^2+(plg{1}.pf.y-plg{1}.pl.y)^2);
                c= sqrt((plg{end}.pf.x-plg{1}.pl.x)^2+(plg{end}.pf.y-plg{1}.pl.y)^2);
                cosTheta= (c^2-b^2-a^2)./(-2*a*b);
                theta=rad2deg(acos(cosTheta));
                if (plg{1}.dis>minDis) || (plg{end}.dis>minDis)
                    if theta>maxTheta
                        plg{1}.pf=plg{end}.pf;%change the first point
                        plg{1}.dis= sqrt( (plg{1}.pl.x-plg{1}.pf.x)^2 + (plg{1}.pl.y-plg{1}.pf.y)^2 );
                        plg{1}.v1=plg{1}.pl.x-plg{1}.pf.x;
                        plg{1}.v2=plg{1}.pl.y-plg{1}.pf.y;
                        plg(end)=[];
                        limit=limit-1;
                    end
                end
            end
        end
        %Getting the angles between segments
        for n=1:size(plg,2)-1
            a= sqrt((plg{n}.pf.x-plg{n}.pl.x)^2+(plg{n}.pf.y-plg{n}.pl.y)^2);
            b= sqrt((plg{n+1}.pf.x-plg{n+1}.pl.x)^2+(plg{n+1}.pf.y-plg{n+1}.pl.y)^2);
            c= sqrt((plg{n}.pf.x-plg{n+1}.pl.x)^2+(plg{n}.pf.y-plg{n+1}.pl.y)^2);
            cosTheta= (c^2-b^2-a^2)./(-2*a*b);
            plg{n}.thetaAS=rad2deg(acos(cosTheta));
        end
        a= sqrt((plg{end}.pf.x-plg{end}.pl.x)^2+(plg{end}.pf.y-plg{end}.pl.y)^2);
        b= sqrt((plg{1}.pf.x-plg{1}.pl.x)^2+(plg{1}.pf.y-plg{1}.pl.y)^2);
        c= sqrt((plg{end}.pf.x-plg{1}.pl.x)^2+(plg{end}.pf.y-plg{1}.pl.y)^2);
        cosTheta= (c^2-b^2-a^2)./(-2*a*b);
        plg{end}.thetaAS=rad2deg(acos(cosTheta));
        %Now plg have the filtered polygon
        cfn=0;%Get the turning angle function of the polygon
        thetaS=[]; thetaAS=[]; disS=[]; taf=[]; cf=[]; cfP=[];
        rthetaS=[]; rthetaAS=[]; rdisS=[]; rtaf=[]; rcf=[]; rcfP=[];
        srtaf=[]; srdisS=[]; sptaf=[]; spdisS=[]; 
        pcf=[]; pcfP=[]; ptaf=[]; pdisS=[];
                        
        maxy1=BWy+1;
        for n=1:size(plg,2)
            theta=rad2deg(atan2(( (maxy1-plg{n}.pl.y) - (maxy1-plg{n}.pf.y) ),(plg{n}.pl.x-plg{n}.pf.x)));
            thetaS(n)=theta; %Angle of the segment (in deg format)
            thetaAS(n)=plg{n}.thetaAS; %Angle between adjacent segments (in deg format)
            disS(n)=plg{n}.dis; %Lenght of the segment
        end
        
        taf(1)=thetaS(1)-thetaS(n);
        for n=2:size(plg,2)
            taf(n)=thetaS(n)-thetaS(n-1);
        end              
        %set the turn angle in range [-180,180]
        for n=1:size(plg,2)
           if(taf(n)>180) 
                taf(n)=taf(n)-360;
               else if (taf(n)<-180)                
                      taf(n)=360+taf(n);
               end
           end
        end
        %sorting the TAF to start in the beginning of a new spicule
        fp=0;
        for n=1:size(taf,2)-1
            if (taf(n)>=0 && taf(n+1)<0)
                fp=n+1;
                initTheta=n;
                break;
            end
        end
        if ( (taf(n+1)>=0) && (taf(1)<0) || (fp==0) )
            fp=1;
            initTheta=n+1;
        end
        if fp==1
            phi=thetaS(fp);%first angle
            rtaf=taf; 
            rdisS=disS;
            rthetaS=thetaS;
            rthetaAS=thetaAS;
        else
            phi=thetaS(fp);
            nmod=0;
            for n=fp:size(taf,2)
                nmod=nmod+1;
                rtaf(nmod)=taf(n);
                rdisS(nmod)=disS(n);
                rthetaS(nmod)=thetaS(n);
                rthetaAS(nmod)=thetaAS(n);
            end
            for n=1:fp-1
                nmod=nmod+1;
                rtaf(nmod)=taf(n);
                rdisS(nmod)=disS(n);
                rthetaS(nmod)=thetaS(n);
                rthetaAS(nmod)=thetaAS(n);
            end
        end
        %constructing the filtered TAF function to be plotted
        rcfn=0;

        for n=1:size(rtaf,2)
            if n==1%rthetaS(fp) means the initial phi in the book
                rcf(1)=thetaS(initTheta); %initTheta is unsorted, then we must use the taf and NOT the sorted rtaf
                ptaf(1)=thetaS(initTheta);%taf to be plotted
                pdisS(1)=disS(initTheta);
                for m=1:round(disS(initTheta))
                 rcfn=rcfn+1;
                 rcfP(rcfn)=rcf(n);
                end
            else
                rcf(n)=rtaf(n-1)+rcf(n-1);
                ptaf(n)=rtaf(n-1);
                pdisS(n)=rdisS(n-1);
                for m=1:round(rdisS(n-1))
                  rcfn=rcfn+1;
                  rcfP(rcfn)=rcf(n);
                end
            end
                          
        end
           rcf(n+1)=taf(initTheta)+rcf(n); %initTheta is unsorted, then we must use the taf and NOT the sorted rtaf
           ptaf(n+1)=taf(initTheta);
           pdisS(n+1)=disS(initTheta);
           rcfn=rcfn+1;
           rcfP(rcfn)=rcf(n+1);%filtered TAF to be plotted
        %Get a new spicules each time that we start with a new decrease in the
        %angle of the segment (taf) differences
        spn=1; %spicules number
        spl=[]; spl(1)=spn;
        for n=2:size(rtaf,2)            
            if rtaf(n-1)>=0 && rtaf(n)<0
                spn=spn+1;
            end
            spl(n)=spn;
        end
        %compute descriptors from filtered TAF
        xTAFn=xTAFn+1;
        [xF,featsF] = getDescriptorsTAF ( rthetaAS, rtaf, rdisS, spl);
        featsfTAF=horzcat(featsfTAF,strcat('taf_maxTheta',num2str(setOfTheta(sot)),'_minDis',num2str(setOfDistances(sod)),'_',featsF));
        xfTAF=horzcat(xfTAF,xF);                       
        %--------------------
        %Getting the compression rate
         Cv=[];
         Cv(1,1)=plg{1}.pf.y; Cv(1,2)=plg{1}.pf.x;%first point of the segment
         for n=1:size(plg,2) %Get the points of the polygon
             Cv(n+1,1)=plg{n}.pl.y; Cv(n+1,2)=plg{n}.pl.x;%last point of the segment
         end
         BWCv = poly2mask(Cv(:,2), Cv(:,1), BWy, BWx);
                           
         %Crp=number of points of the new polygon / number of points of the original poligon
          Crp=(size(Cv,1)-1)./(size(Co,1)-1); %compressi??n ratio        
          [hausdd] = getHausddDistance(BWCo, BWCv);
          interCoCv=sum(sum((BWCo&BWCv)));
          joinCoCv=sum(sum((BWCo|BWCv)));
          JI=interCoCv/joinCoCv;%jaccard index
          JD=(joinCoCv-interCoCv)/joinCoCv;%jaccard distance -> is equal to JD=1-JI
          %Average radial distance and m??ximum radial distance
          [aveRD, maxRD] = getNRadialDistance(BWCo, BWCv);
          xPV=[Crp, hausdd, JI, aveRD, maxRD];
          featsPV={'Cpr', 'hausdd', 'JI', 'JD', 'aveRD', 'maxRD'};
          featsPM=horzcat(featsPM,strcat('taf_maxTheta',num2str(setOfTheta(sot)),'_minDis',num2str(setOfDistances(sod)),'_',featsPV));
          xPM=horzcat(xPM,xPV,JI);  
        %--------------------
        %construct smoothed TAF grouping the elements with the same sign
        nmod=0;
        m=1;
        while m<=size(rtaf,2)
            nmod=nmod+1;
            if m==size(rtaf,2)                
                srtaf(nmod)=rtaf(m);
                srdisS(nmod)=rdisS(m);
                m=m+1;
            else                           
            n=m+1;
            while (n<=size(rtaf,2))&&(sign(rtaf(m))==sign(rtaf(n)))
             n=n+1;   
            end                                 
            srtaf(nmod)=mean(rtaf(m:n-1));
            srdisS(nmod)=mean(rdisS(m:n-1));
            m=n;%This is the next element with different sign 
            end
        end
        %identify spiculations in the smoothed TAF
        %the regions joined by an increased angle
        spn=1; %spicules number
        spl=[]; spl(1)=spn;
        for n=2:size(srtaf,2)            
            if srtaf(n-1)>=0 && srtaf(n)<0
                spn=spn+1;
            end
            spl(n)=spn;
        end
        %identify the regions joined by a decreased angle
        
        spn=1; %spicules number
        spl2=[]; spl2(1)=1;
        for n=2:size(srtaf,2)            
            if srtaf(n-1)<0 && srtaf(n)>0
                spn=spn+1;
            end
            spl2(n)=spn;
        end
          if srtaf(end)>0 && srtaf(1)<0
            spl2(end)=1;    
          end
        
        [xS,featsS] = getDescriptorsSTAF ( srtaf, srdisS, spl, spl2);
        featsSTAF=horzcat(featsSTAF,strcat('taf_maxTheta',num2str(setOfTheta(sot)),'_minDis',num2str(setOfDistances(sod)),'_',featsS));
        xSTAF=horzcat(xSTAF,xS); 
        %--------------------
        %construct smoothed TAF to be plotted grouping the elements with the same sign
        nmod=0;
        m=1; nmod=nmod+1; %ptaf(1) is the phi angle, the orientation angle of the polygon
        sptaf(nmod)=ptaf(m);
        spdisS(nmod)=pdisS(m);
        m=2;
        while m<=size(ptaf,2)
            nmod=nmod+1;
            if m==size(ptaf,2)                
                sptaf(nmod)=ptaf(m);
                spdisS(nmod)=pdisS(m);
                m=m+1;
            else                           
            n=m+1;
            while (n<=size(ptaf,2))&&(sign(ptaf(m))==sign(ptaf(n)))
             n=n+1;   
            end                                 
            sptaf(nmod)=mean(ptaf(m:n-1));
            spdisS(nmod)=mean(pdisS(m:n-1));
            m=n;%This is the next element with different sign 
            end
        end
        rcfn=0;
        for n=1:size(sptaf,2)
           pcf(n)=sptaf(n);                                 
            for m=1:round(spdisS(n))
              rcfn=rcfn+1;
              pcfP(rcfn)=pcf(n);%smoothed TAF to be plotted
            end                        
        end
        %--------------------
    end
end

if numel(idxStats) == 3
    varargout{1} = xfTAF;
    varargout{2} = featsfTAF;
    varargout{3} = xSTAF;
    varargout{4} = featsSTAF;
    varargout{5} = xPM;
    varargout{6} = featsPM;    
elseif (numel(idxStats) == 1) && (idxStats == 1)
    varargout{1} = xfTAF;
    varargout{2} = featsfTAF;
elseif (numel(idxStats) == 1) && (idxStats == 2)
    varargout{1} = xSTAF;
    varargout{2} = featsSTAF;
elseif (numel(idxStats) == 1) && (idxStats == 3)
    varargout{1} = xPM;
    varargout{2} = featsPM;
elseif (numel(idxStats) == 2) && (idxStats(1) == 1) && (idxStats(2) == 2)
    varargout{1} = xfTAF;
    varargout{2} = featsfTAF;
    varargout{3} = xSTAF;
    varargout{4} = featsSTAF;    
elseif (numel(idxStats) == 2) && (idxStats(1) == 1) && (idxStats(2) == 3)
    varargout{1} = xfTAF;
    varargout{2} = featsfTAF;
    varargout{3} = xPM;
    varargout{4} = featsPM;  
elseif (numel(idxStats) == 2) && (idxStats(1) == 2) && (idxStats(2) == 3)
    varargout{1} = xSTAF;
    varargout{2} = featsSTAF;
    varargout{3} = xPM;
    varargout{4} = featsPM;   
end

%***********************************************************************
function [disHaussdd] = getHausddDistance(BWL, BWP)
%boundary of the lesion
junk = bwboundaries(BWL);
cBW  = junk{1};
yBW  = cBW(:,1); xBW = cBW(:,2);
%boundary of the 
junk = bwboundaries(BWP);
cP = junk{1};
yP  = cP(:,1); xP = cP(:,2);

D1=dist(cP,cBW');
[Dmp,ind] = min(D1,[],2);
%distancia del poligono a la lesi??n
[disHaussdd, ind2]=max(Dmp);


%Get the radial distance (normalized by the size of the lesion) between the lesion boundary and the polygonal
%model
% BWL - boundary of the lesion
% BWP - boundary of the polygon
%***********************************************************************
function [disn, maxdis] = getNRadialDistance(BWL, BWP)
% Extrae y parametriza contorno del tumor 
junk = bwboundaries(BWL);
cBW  = junk{1};
yBW  = cBW(:,1); xBW = cBW(:,2);
props = regionprops(BWL,'Centroid');
xc = props.Centroid(1);
yc = props.Centroid(2);
% Vectores unitarios de BW de la lesi??n
rBW = [xBW-xc yBW-yc];
nBW = sqrt(sum(rBW.^2,2));
uBW = rBW./(repmat(nBW,1,2)+eps);
thetaBW = atan2(rBW(:,2),rBW(:,1));
% Extrae coordenadas del pol??gono 
junk = bwboundaries(BWP);
cP = junk{1};
yP  = cP(:,1); xP = cP(:,2);
% Empareja coordenadas del pol??gono
% % D1 = dist([xBW yBW],[xP yP]');
% % [~,ind] = min(D1,[],2);
% % xP = xBW(ind);
% % yP = yBW(ind);
% % 
% % CP = roipoly(BWL,xP,yP);
% % junk = bwboundaries(CP);
% % cP = junk{1};
% % yP  = cP(:,1); xP = cP(:,2);
% Vectores unitarios de BWP del pol??gono
rP = [xP-xc yP-yc];
nP = sqrt(sum(rP.^2,2));
uP = rP./(repmat(nP,1,2)+eps);
thetaP = atan2(rP(:,2),rP(:,1));
% Distancia entre vectores unitarios
D = dist(uP,uBW');
[~,ind] = min(D,[],2);
% Correspondencia entre puntos de BW y puntos en CH con la orientacion
% mas proxima
mdBW = cBW(ind,:);
% Distancia entre puntos del contorno de la lesi??n y los puntos del
% polygono con la orientaci??n m??s pr??xima
%%---------------------------aqui voy
dis = diag(dist(cP,mdBW'));
maxdis=max(dis);
disn=sum(dis)/(size(cBW,1));
%Spiculation index (0 <= SItf <= 1) -SItf is higher for malignant lesions and lower for benign lesions
%thetaS= angle of the segment, disS=lenght of the segment,
%spl=identificator index of segments
%thetaAS=angle between adjacent segments
%***********************************************************************
function [x,feats] = getDescriptorsTAF ( thetaAS, rtaf, rdisS, spl)
%%Getting the spicules identificators
uspl=unique(spl);
n=0;
for m=1:size(uspl,2)
    %each spiculation at least must have at leat one decrease and one increase in the taf
    if (sum(spl==uspl(m)&rtaf>0))>0 && (sum(spl==uspl(m)&rtaf<0))>0
        n=n+1;
        
        %remove the last segment, because the angles between segments are taken with the next segments
        %so in this way we only keep the angles inside the segments of the spiculation
        sE=spl==uspl(m); 
        for x=size(sE,2):-1:1
            if(sE(x)==1)
                sE(x)=0;
                break;
            end
        end
        lp(n)=sum(rdisS(spl==uspl(m)));%sum of distances of all sections
        mptheta(n)=mean(thetaAS(sE)); %mean of internal angles of the sections of the spiculation
        ptheta(n)=mean(thetaAS(sE&thetaAS<=mptheta(n))); %internal angle of the spiculation
        
    end
end
if n>0
    SItf=sum((1+cos(deg2rad(ptheta))).*lp)./(2.*sum(lp));
else
    SItf=0;
end
x=[SItf];
feats={'SItf'};

%CXta - Index of convexity from the smoothed taf, CXta is in the range
%[0,1], where a convex contour is equal to 1
%VRta - index of concave regions, range [0,1], for convex contour VRta is 0
%XRta - index of convex regions, range [0,1], for convex contour XRta is 1
%SIta - spiculation index from smoothed taf
%***********************************************************************
function [x,feats] = getDescriptorsSTAF ( srtaf, srdisS, spl, spl2)
%%Getting the spicules identificators
uspl=unique(spl);
n=0;
for m=1:size(uspl,2)
    %each spiculation at least must have at leat one decrease and one increase in the taf
    %and must have only two sections
    if (sum(spl==uspl(m)&srtaf>0))>0 && (sum(spl==uspl(m)&srtaf<0))>0 && sum(spl==uspl(m))==2
        n=n+1;
        sE=spl==uspl(m); %searching for the index of the last section, which also is the positive section
        for x=size(sE,2):-1:1 %x will be the positive section
            if(sE(x)==1)                
                break;
            end
        end
        %SIta features
        %thetaP=180-|Tc(sp+1)-Tc(sp)|
        %given that (Tc(sp+1)-Tc(sp))=is the positive turning angle of sp+1 or
        %srtaf(sp+1), where srtaf(sp+1)>=0 is the positive section
        %then ptheta=180-abs(srtaf(sp+1))
        ptheta(n)=180-abs(srtaf(x));  %this is the internal angle of the spicule n
        lp(n)=sum(srdisS(spl==uspl(m))); %sum of the lenghts of the sides of the spicule n
        %XRta features
        %This pair of sections are joined by an increasing angle
        XRtheta(n)=180-abs(srtaf(x)); 
        XRlp(n)=sum(srdisS(spl==uspl(m)));                         
    end
end
%Getting the regions joined by a drop angle
uspl2=unique(spl2);
n2=0;
for m=1:size(uspl2,2)
    if (sum(spl2==uspl2(m)&srtaf>0))>0 && (sum(spl2==uspl2(m)&srtaf<0))>0 && sum(spl2==uspl2(m))==2
        n2=n2+1;
        sE2=spl2==uspl2(m); %searching for the index of the last section, which also is the negative section
        for x2=size(sE2,2)-1:-1:1 %x will be the negative section, the search must start in size(sE2,2)-1 because the position size(sE2,2) is positive
            if(sE2(x2)==1)                
                break;
            end
        end
        if (x2==0)&&sE2(size(sE2,2)==1)
            x2=size(sE2,2);
        end
        %VRta features
        %This pair of sections are joined by a decreasing angle
        VRtheta(n2)=180-abs(srtaf(x2)); 
        VRlp(n2)=sum(srdisS(spl2==uspl2(m))); 
    end
end

if n>0
    SIta=(sum((1+cos(deg2rad(ptheta))).*lp)) ./ (2*sum(lp));
    XRta=1-(sum((1+cos(deg2rad(XRtheta))).*XRlp)) ./ (2*sum(XRlp));
else
    SIta=0;
    XRta=1;
end

if n2>0
    VRta=(sum((1+cos(deg2rad(VRtheta))).*VRlp)) ./ (2*sum(VRlp));
else
    VRta=0;
end

CXta=(XRta/2)+((1-VRta)/2);

x=[SIta VRta XRta CXta];
feats={'SIta' 'VRta' 'XRta' 'CXta'};
%***********************************************************************
function r = max_inner_circle(x,y)
warning off;
% make a voronoi diagram
[vx,vy]=voronoi(x,y);
% find voronoi nodes inside the polygon [x,y]
Vx=vx(:);
Vy=vy(:);
% Here, you could also use a faster version of inpolygon
IN=inpolygon(Vx,Vy, x,y);
ind=find(IN==1);
Vx=Vx(ind);
Vy=Vy(ind);
% maximize the distance of each voronoi node to the closest node on the
% polygon.
minDist=0;
for i=1:length(Vx)
    dx=(Vx(i)-x);
    dy=(Vy(i)-y);
    r=min(dx.*dx+dy.*dy);
    if (r>minDist)
        minDist=r;
    end
end
r = sqrt(minDist);