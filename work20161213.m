%% load in movie of motile microtubules
obj = VideoReader('tubgfp2_small.mov')
outfilename = 'tubgfp2_small.tiff'
for img = 1:obj.NumberOfFrames;
    %filename = strcat('frame',num2str(img),'.jpg');
    b = read(obj,img);
    imshow(b)
    drawnow
    imwrite(b,outfilename,'WriteMode','append','Compression','none');
end