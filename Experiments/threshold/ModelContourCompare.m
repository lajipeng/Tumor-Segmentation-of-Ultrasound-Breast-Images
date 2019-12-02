function ModelContourCompare(subset)
%-- Author：Wang, Peng. Fudan University 
%-- Date：2019/11/30
%-- Fuction：show the contour of segmentation result from diffirent model
%% loading and initializating image
load('C:/Github/医学成像技术/project2/Tumor-Segmentation-of-Ultrasound-Breast-Images/Input.mat')
load('C:/Github/医学成像技术/project2/Tumor-Segmentation-of-Ultrasound-Breast-Images/GT.mat')
%% imshow 
original_image = Input{subset};
GroundTrue = logical(imfill(GT{subset},'holes'));
% [A autosegment_image] = autosegment(original_image)
horsh_image = horsch(original_image,GroundTrue);
wtsdsegment_image = wtsdsegment(original_image,GroundTrue);
threshold_image = thresholdsegment(original_image);
%%
[B1,L1] = bwboundaries(GroundTrue,'noholes');
[B2,L2] = bwboundaries(threshold_image,'noholes');
[B3,L3] = bwboundaries(horsh_image,'noholes');
[B4,L4] = bwboundaries(autosegment_image,'noholes');
[B5,L5] = bwboundaries(wtsdsegment_image,'noholes');

boundary1 = B1{1}; 
boundary2 = B2{1}; 
boundary3 = B3{1}; 
boundary4 = B4{1};
boundary5 = B5{1}; 
%%
figure;
subplot 151; imshow(original_image,'InitialMagnification', 'fit'); xlabel('(a)');
subplot 152; imshow(original_image,'InitialMagnification', 'fit'); xlabel('(b)');
hold on; 
plot(boundary1(:, 2), boundary1(:, 1), 'r', 'LineWidth', 2); 
plot(boundary2(:, 2), boundary2(:, 1), 'g', 'LineWidth', 2);
subplot 153; imshow(original_image,'InitialMagnification', 'fit'); xlabel('(c)');
hold on; 
plot(boundary1(:, 2), boundary1(:, 1), 'r', 'LineWidth', 2); 
plot(boundary3(:, 2), boundary3(:, 1), 'g', 'LineWidth', 2);
subplot 154; imshow(original_image,'InitialMagnification', 'fit'); xlabel('(d)');
hold on; 
plot(boundary1(:, 2), boundary1(:, 1), 'r', 'LineWidth', 2); 
plot(boundary4(:, 2), boundary4(:, 1), 'g', 'LineWidth', 2);
subplot 155; imshow(original_image,'InitialMagnification', 'fit'); xlabel('(e)');
hold on; 
plot(boundary1(:, 2), boundary1(:, 1), 'r', 'LineWidth', 2); 
plot(boundary5(:, 2), boundary5(:, 1), 'g', 'LineWidth', 2);
end