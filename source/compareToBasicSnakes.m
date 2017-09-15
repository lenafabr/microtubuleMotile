%% load in movie of motile microtubules
obj = VideoReader('tubA_GFP_2mic_linear_cap15.mp4')
%%
img = read(obj,17);
%imshow(img,[],'InitialMagnification','fit');
img = double(img(:,:,2));
%imshow(img,[],'InitialMagnification','fit');
img(:) = (img(:)-min(img(:)))/(max(img(:))-min(img(:)));
%% adjust image brightness
img2 = imgaussfilt(img,2,'FilterSize',3);
img3 = imadjust(img2,[0.1,0.6],[0,1],0.7);
imshow(img3,[],'InitialMagnification','fit');
%% crop image
croprect = [563 346  381  366];
img3 = imcrop(img3,croprect);
%%
imshow(img3,[],'InitialMagnification','fit');


%% using Kloon code "BasicSnakes"
tic
Options = struct();
Options.Wline = 1;
Options.Wedge = 0;
Options.Wterm = 0;
Options.Sigma1 = 2;
Options.Sigma2 = 2;
% get gradient of energy function
Eext = ExternalForceImage2D(img3,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1);

Fx=ImageDerivatives2D(Eext,Options.Sigma2,'x');
Fy=ImageDerivatives2D(Eext,Options.Sigma2,'y');
Fext(:,:,1)=Fx*2*Options.Sigma2^2;
Fext(:,:,2)=Fy*2*Options.Sigma2^2;

%% plot gradient field
[X,Y] = meshgrid(1:size(img3,2),1:size(img3,1));
%F1 = Fext(:,:,1); F2 = Fext(:,:,2);
imshow(Eext)
hold all
quiver(X,Y,Fext(:,:,2),Fext(:,:,1),2)
hold off

%% Gradient vector flow
Mu = 0.2;
Iterations = 100;
Sigma = Options.Sigma2;
FextGVF=GVFOptimizeImageForces2D(Fext, Mu, Iterations, Sigma);
toc
%% plot gradientGVFeshgrid(1:size(img3,2),1:size(img3,1));
%F1 = Fext(:,:,1); F2 = Fext(:,:,2);
imshow(img3)
hold all
quiver(X,Y,FextGVF(:,:,2),FextGVF(:,:,1),4)
%quiver(X,Y,Fext(:,:,2),Fext(:,:,1),4)
hold off

%%
P = [247          180
          194           93
          237           71];
hold all; plot(P(:,1),P(:,2),'.-')
      
%%
% Get image force on the contour points
kappa = 1;
Fext1(:,1)=kappa*interp2(FextGVF(:,:,1),P(:,1),P(:,2));
Fext1(:,2)=kappa*interp2(FextGVF(:,:,2),P(:,1),P(:,2));
% Interp2, can give nan's if contour close to border
Fext1(isnan(Fext1))=0;

hold all
quiver(P(:,1),P(:,2),Fext1(:,2),Fext1(:,1),1)

