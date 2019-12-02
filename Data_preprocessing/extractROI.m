%第一步：从图片中选取矩形框区域
I = imread('4.jpg');
[A,rect] = imcrop(I);
imshow(A);
rect
%第二步：根据rect确定：在原图中绘制的矩形的坐标，注意rect的格式[m n l k]->[(m,n) (m+l,n+k)]->[(n,m) (n+k,m+l)]    

 

I2 = imcrop(I,rect);
figure(2);imshow(I2);  
                                    
                   