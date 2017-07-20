img = imread('frame1.png');
img = double(img(:,:,2));
imshow(img,[],'InitialMagnification','fit');

%% adjust image brightness
%img2 = imgaussfilt(img,2,'FilterSize',3);
img3 = imadjust(img,[0.1,0.6],[0,1],0.7);
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
     [u,v] = GVF(img3, 0.2, 40); 
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