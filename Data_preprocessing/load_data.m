clc;
clear;
%% this to read avi by using videoread to get every frame
video = VideoReader('C:\Github\医学成像技术\project2\project_2乳腺肿瘤.avi');
nFrames = video.NumberOfFrames;   %get frames
H = video.Height;     %get Hight
W = video.Width;      %get Width
Rate = video.FrameRate;
Ultra_images = {};
% Preallocate movie structure.
mov(1:nFrames) = struct('cdata',zeros(H,W,3,'uint8'),'colormap',[]);
%read one frame every time
for i = 1:nFrames
    mov(i).cdata = read(video,i);
    P = mov(i).cdata;
    disp('当前播帧数：'),disp(i);
    figure(1),
    imshow(P),title('原始图片');
    Ultra_images{i} = P;
end
%% this to read ground true 
load('C:\Github\医学成像技术\project2\Tumor-Segmentation-of-Ultrasound-Breast-Images\src\project_2_manual segmentation')
for i = 1:nFrames
    disp('当前播帧数：'),disp(i);
    figure(2),
    imshow(BW(:,:,i),[]),title('原始图片');
end
%% this to show groud true in original images
load('C:\Github\医学成像技术\project2\Tumor-Segmentation-of-Ultrasound-Breast-Images\src\project_2_manual segmentation.mat')
load('C:\Github\医学成像技术\project2\Tumor-Segmentation-of-Ultrasound-Breast-Images\src\Project_2_OriginImages.mat')
for i = 1:nFrames
    disp('当前播帧数：'),disp(i);
    figure(3),
    imshowpair(Ultra_images{i},BW(:,:,i),'diff'),title('原始图片');
end
%% this is to extract the region of interst
load('C:\Github\医学成像技术\project2\Tumor-Segmentation-of-Ultrasound-Breast-Images\src\project_2_manual segmentation.mat')
load('C:\Github\医学成像技术\project2\Tumor-Segmentation-of-Ultrasound-Breast-Images\src\Project_2_OriginImages.mat')
rect = [359,162,300,220];
Input = {};
GT  = {};
for i = 1:nFrames
    ground_true = imcrop(BW(:,:,i),rect);
    origin_image = imcrop(Ultra_images{i},rect);
    disp('当前播帧数：'),disp(i);
    figure(4),
    imshowpair(origin_image,ground_true,'diff'),title('Region of Interest');
    GT{i} = ground_true;
    Input{i}= origin_image(:,:,1);
    filename = sprintf('Input%d.bmp',i);
    imwrite(filename,origin_image(:,:,1));
end
