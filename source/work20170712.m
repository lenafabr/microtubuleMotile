
%% load in movie of motile microtubules
obj = VideoReader('tubA_GFP_2mic_linear_cap15.mp4')
outfilename = 'movie15.tiff'
%%
for img = 1:obj.NumberOfFrames;
    %filename = strcat('frame',num2str(img),'.jpg');
    b = read(obj,img);
    imshow(b)
    drawnow
    imwrite(b(:,:,2),outfilename,'WriteMode','append','Compression','none');
end
%%
img = read(obj,17); img = img(:,:,2);
imshow(img,[],'InitialMagnification','fit');

%% adjust image brightness
img2 = imgaussfilt(img,2,'FilterSize',3);
img3 = imadjust(img2,[0.1,0.6],[0,1],0.7);
imshow(img3,[],'InitialMagnification','fit');
%% crop image
croprect = [563 346  381  366];
img3 = imcrop(img3,croprect);
%%
imshow(img3,[],'InitialMagnification','fit');

%% pick a starting pixel
pstart = [167,253];
hold all
plot(pstart(1), pstart(2),'y.')
hold off
%% localized radon transform at the starting pixel
[X,Y] = meshgrid(1:size(img3,1),1:size(img3,2));
img3d = double(img3);

L=100; % length to look over
theta = 22*pi/180;
dvec = [cos(theta),sin(theta)];% direction
lams = -L/2:1:L/2;
linepts = [dvec(1)*lams'+pstart(1),dvec(2)*lams'+pstart(2)];
hold all
plot(linepts(:,1),linepts(:,2),'c')
hold off
%%
L=100; % length to look over
lams = -L/2:1:L/2;

ntheta = 180;
thetavals = linspace(0,179,ntheta);
Rvals = zeros(size(thetavals));
for tc = 1:length(thetavals)
    theta = thetavals(tc)*pi/180;
    dvec = [cos(theta),sin(theta)];% direction
    linepts = [dvec(1)*lams'+pstart(1),dvec(2)*lams'+pstart(2)];
    Rvals(tc) = sum(interp2(img3d,linepts(:,1),linepts(:,2)));
end
plot(thetavals,Rvals)

%% test radon transform function

[maxR,maxth,Rvals,thetavals] = maxFiniteRadon(img3d,pstart,L);