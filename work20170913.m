opt = struct();
opt.Lradon = 20; % length to use for finite radon transform
opt.Ltrace = 10; % separation between points for the tracing
% cutoff for radon transform magnitude to call filament end
% relative to starting point magnitude
% cutoff in theta change to call filament end
% in degrees
opt.relthcutoff = 60;

% max points to trace in each direction
opt.maxpts = 100;

% display tracing during calculation
opt.dodisplay = 1;

opt.thEnergy = 100;
opt.thRelEnergy = 0.6;
opt.thAngle = 40;% was 70
opt.thAngleDiff = 25;
opt.nPeaks = 4;
opt.path_residue = 16; % path closeness th r^2) 

obj = VideoReader('tubA_GFP_2mic_linear_cap15.mp4')
outfilename = 'movie15.tiff'

%%
img = read(obj,17); img = img(:,:,2);
imshow(img,[],'InitialMagnification','fit');

% adjust image brightness
img2 = imgaussfilt(img,2,'FilterSize',3);
img3 = imadjust(img2,[0.1,0.6],[0,1],0.7);
imshow(img3,[],'InitialMagnification','fit');
% crop image
croprect = [563 346  381  366];
img3 = imcrop(img3,croprect);
img3d = double(img3);
%% sat the images to reduce dynamic range since intensity at core plays little role.
% img3d(find(img3d>150))=150;

%
imshow(img3d,[],'InitialMagnification','fit');


% pick a starting pixel
% pstart = [167,253];
% pstart = [204,110];
% pstart = [231,228];
% pstart = [204,81];
% pstart = [198,101];
% pstart = [276,92];
% pstart = [136, 249];
 pstart = ginput(1);
% pstart = [248, 120];
% pstart = [72, 248];
% pstart = [227,71];
% pstart = [275.659,92.5800];

hold all
plot(pstart(1), pstart(2),'y.','MarkerSize',24)
hold off

position = pstart;
init_Intensity = 5000;
min = zeros;
angle2compare = 0;
global path_vertexs;
global max_pos;
global max_compare;
path_vertexs = pstart;
max_pos = zeros;
max_compare = zeros; 

RadonTree(img3d, pstart, init_Intensity, position, angle2compare, opt);

disp(max_pos);
hold all;
plot(max_pos(:,1),max_pos(:,2),'c', 'lineWidth', 5);
hold off;
res = 22050;
len = 0.5 * res;

disp("done")



