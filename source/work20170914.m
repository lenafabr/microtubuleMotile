obj = VideoReader('../tubA_GFP_2mic_linear_cap15.mp4')

%% load in image
img = read(obj,26); 
%imshow(img,[],'InitialMagnification','fit');

%img = imread('frame1.png');
img = double(img(:,:,2));
%imshow(img,[],'InitialMagnification','fit');
img(:) = (img(:)-min(img(:)))/(max(img(:))-min(img(:)));
imshow(img,[],'InitialMagnification','fit');
% %% adjust image brightness
% img2 = imgaussfilt(img,2,'FilterSize',3);
% img3 = imadjust(img2,[0.1,0.6],[0,1],0.7);
% imshow(img3,[],'InitialMagnification','fit');
%% crop image
croprect = [563 346  381  366];
img3 = imcrop(img,croprect);
imshow(img3,[])

%% calculate GVF field
[FextGVF,Fext] = getGVFfield(img3);

% display gvf field
[X,Y] = meshgrid(1:size(img3,2),1:size(img3,1));
%F1 = Fext(:,:,1); F2 = Fext(:,:,2);
imshow(Eext,[])
hold all
quiver(X,Y,FextGVF(:,:,2),FextGVF(:,:,1),2)
hold off

%% set initial guess
ptinit = [221.91       234.03
       177.13       186.49
       183.76       159.95
       196.48       140.05];


%% set up chain object
param = arclenparam(ptinit');

% constraints
% fix end positions
chain.pos0 = ptinit(1,:)'; 
chain.posf = ptinit(end,:)';
% do not fix end tangents
chain.fixtan0 = 0; chain.fixpos0 = 0; 
chain.fixtanf = 0; chain.fixposf = 0;
chain.tan0 = [1,0];
chain.tanf = [1,0];


chain.nseg = 20;
chain.nbead = chain.nseg-1;
chain.ncrd = 2*chain.nbead;

chain.len = param(end); % overall length of chain
if (chain.fixpos0)
    chain.ls = chain.len/chain.nseg; % ground state segment length
else
    chain.ls = chain.len/(chain.nseg-2);
end

% energetics
chain.lp = 1; % persistence length (nm)
chain.lstretch = 1; % stretch modulus (kT/nm)


   
%% try running filament evolution over time
savecoords = evolveFilament(obj,chain,ptinit,struct('frames',26:1:100,'croprect',croprect,...
    'stepsize',0.5,'steptol',1e-4,'extscl',1,'getgvf',0))