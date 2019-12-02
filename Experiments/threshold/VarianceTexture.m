function [textureIm glcm] = VarianceTexture(im)

imGray = mat2gray(im);
[ren col] = size(imGray);

textureIm = zeros(ren,col);

r = 6;     % Adjust for desired window size

M = ren+2*r;
N = col+2*r;

imPad = padarray(imGray,[r r],'symmetric');

offsets = [[0 -1 -1 -1]' , [1 1 -1 -1]'];

for n = 1+r:N-r
    for m = 1+r:M-r
        
        window = imPad(m+(-r:r),n+(-r:r));

        glcm = graycomatrix(window,'GrayLimits',[0 1],'NumLevels',16,'Offset',offsets);
        glcm = (glcm(:,:,1) + glcm(:,:,2) + glcm(:,:,3) + glcm(:,:,4))/4;

        variance = var(glcm(:));

        textureIm(m-r,n-r) = variance;
        
    end
end

textureIm = mat2gray(textureIm);
% figure, imshow(textureIm8), title('8 numlevels');
