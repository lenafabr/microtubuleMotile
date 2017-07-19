%% load in movie of motile microtubules
obj = VideoReader('tubA_GFP_2mic_linear_cap15.mp4')
outfilename = 'movie15.tiff'
for img = 1:obj.NumberOfFrames;
    %filename = strcat('frame',num2str(img),'.jpg');
    b = read(obj,img);
    imshow(b)
    drawnow
    imwrite(b(:,:,2),outfilename,'WriteMode','append','Compression','none');
end
%% read images from tiff file
img = imread('movie15.tiff',13);
imshow(img,[],'InitialMagnification','fit')
%%
img = read(obj,17); img = img(:,:,2);
imshow(img,[],'InitialMagnification','fit');

%% show 2 images together
fr = 40;
imgA = read(obj,fr); imgA = imgA(:,:,2);
imgB = read(obj,fr+1); imgB = imgB(:,:,2);
imshowpair(imgA,imgB,'falsecolor')
%% Gaussian filter
img2 = imgaussfilt(img,2,'FilterSize',5);
imshow(img2,[],'InitialMagnification','fit');

%% edge detection
BW = edge(img2,'Prewitt');
imshow(BW,'InitialMagnification','fit')

%% adjust image brightness
img2 = imgaussfilt(img,2,'FilterSize',3);
img3 = imadjust(img2,[0.1,0.6],[0,1],0.7);
imshow(img3,[],'InitialMagnification','fit');

%% crop image
croprect = [563 346  381  366];
img3 = imcrop(img3,croprect);
%%
imshow(img3,[]);

%% dilation
SE = strel('disk',2,0);
img4 = imdilate(img3,[SE90,SE0]);
%img4 = imdilate(img3,SE);
imshow(img4,[],'InitialMagnification','fit');


%% thresholding

level = graythresh(img3);
BW = im2bw(img3,level*0.5);
imshow(BW,'InitialMagnification','fit')

%% img open
SE90 = strel('line',5,90);
SE0 = strel('line',5,0);
SE = strel('disk',2,0);
BW2 = imopen(BW,SE);
imshow(BW2,[],'InitialMagnification','fit');

%%
imshowpair(img3,BW2,'falsecolor')


