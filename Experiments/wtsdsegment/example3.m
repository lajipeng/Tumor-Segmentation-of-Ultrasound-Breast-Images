%-- Author：Wang, Peng. Fudan University 
%-- Date：2019/11/30
%-- Fuction：Breast Tumor Segmentation of Ultra Images 
%-- Reference:
%   https://www.springerprofessional.de/en/busat-a-matlab-toolbox-for-breast-ultrasound-image-analysis/
%   https://www.tamps.cinvestav.mx/~wgomez/downloads.html
%% loading and initializating image
load('C:\Github\医学成像技术\project2\Tumor-Segmentation-of-Ultrasound-Breast-Images\Input.mat')
load('C:\Github\医学成像技术\project2\Tumor-Segmentation-of-Ultrasound-Breast-Images\GT.mat')
%% segmentation, compute dice, show the comparasion between the input and ground truth 
nFrames = 95;
wtsdsegment_dice = [];
for i=1:nFrames   
    Smanual = logical(imfill(GT{i},'holes'));
    I = Input{i};
    S = wtsdsegment(I,Smanual);
    wtsdsegment_dice = [wtsdsegment_dice,dice(S,Smanual)];
    imshowpair(Smanual,S,'falsecolor');
    h=findobj(gcf,'type','image');
    result=get(h,'CData');
    filename =  sprintf('results/comparasion/wtsdsegment%d.jpg',i);
    imwrite(result,filename);
end
wtsdsegment_dice_mean = mean(wtsdsegment_dice);
save('results/dice/wtsdsegment_dice.mat','wtsdsegment_dice','wtsdsegment_dice_mean');
%% compute time consuming 
t1=clock;
for i=1:nFrames   
    Smanual = logical(imfill(GT{i},'holes'));
    I = Input{i};
    S = wtsdsegment(I,Smanual);
end
t2=clock;
time = etime(t2,t1)/nFrames;
disp(['etime程序总运行时间：',num2str(time)]);
save('results/time/time.mat','time');