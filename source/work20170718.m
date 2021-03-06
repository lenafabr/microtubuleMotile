img = imread('frame1.png');
img = double(img(:,:,2));
%imshow(img,[],'InitialMagnification','fit');
img(:) = (img(:)-min(img(:)))/(max(img(:))-min(img(:)));
imshow(img,[],'InitialMagnification','fit');
%% adjust image brightness
img2 = imgaussfilt(img,2,'FilterSize',3);
img3 = imadjust(img2,[0.1,0.6],[0,1],0.7);
imshow(img3,[],'InitialMagnification','fit');
%% crop image
croprect = [563 346  381  366];
img3 = imcrop(img,croprect);
%%
imshow(img3,[],'InitialMagnification','fit');

%%

[xinit,yinit] = snakeinit(0.1);

%%
% Compute the GVF of the edge map f
     disp(' Compute GVF ...');
     mu=0.2;
     [u,v] = GVF(img3, mu, 40); 
     disp(' Nomalizing the GVF external force ...');
     mag = sqrt(u.*u+v.*v);
     px = u./(mag+1e-10); py = v./(mag+1e-10); 
     
  %% display the GVF 
 imshow(img3,[]);
 [X,Y] = meshgrid(1:size(img3,2),1:size(img3,1));
 hold all
 s=1;
  quiver(X(1:s:end,1:s:end),Y(1:s:end,1:s:end),px(1:s:end,1:s:end),py(1:s:end,1:s:end));
  axis off; axis equal; axis 'ij';     % fix the axis 
  title('normalized GVF field');
         hold off
         
%% using Kloon code "BasicSnakes"
tic
Options = struct();
Options.Wline = 1;
Options.Wedge = 0;
Options.Wterm = 0;
Options.Sigma1 = 4;
Options.Sigma2 = 4;
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
quiver(X,Y,Fext(:,:,2),Fext(:,:,1),3)
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
quiver(X,Y,Fext(:,:,2),Fext(:,:,1),4)
hold off
