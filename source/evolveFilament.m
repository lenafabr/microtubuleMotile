% evolve filament over time through many frames
% starting with a rough guess of the filament in the first frame
function savecoords = evolveFilament(obj,chain,ptinit,options)
% obj is a videoreader object
% ptinit is an nx2 rough set of initial points describing the filament

% which frames to work on
opt.frames = 1:ceil(obj.duration*obj.framerate);

% which color column to read
opt.col = 2;

% crop rectangle
opt.croprect = [];

% options for optimizing a filament
opt.stepsize = 0.2;
opt.nstep = 10000;
%coords0 = chain.coords;
opt.displayevery = inf;
opt.extscl = 0.5;
opt.getgvf = 1;

% image brightness adjustmant
%opt.adjscl = [0.1,0.6];
%opt.gamma = 0.7;
opt.adjscl = [];
opt.gamma = 1;

% display after optimization is finished
opt.displayafteropt = 1;
opt.displaybeforeopt = 0;

% tolerance in rms step size per bead
% to define convergence
opt.steptol = 1e-4;

if (nargin>3)
    opt = copyStruct(options,opt)
end

% read through initial frames
% for fc = 1:opt.fmin-1
%     obj.readFrame;
% end

%% interpolate initial snake
% parameterize by arc length
param = arclenparam(ptinit');
% interpolate initial snake points
if (chain.fixpos0)
    paramint = linspace(param(1),param(end),chain.nbead+2);
    pts = interp1(param,ptinit,paramint,'linear')
    chain.coords = reshape(pts(2:end-1,:)',2*chain.nbead,1)
else
    paramint = linspace(param(1),param(end),chain.nbead);
    pts = interp1(param,ptinit,paramint,'linear')
    chain.coords = reshape(pts(:,:)',2*chain.nbead,1)
end
%hold all
%plot(pts(:,1),pts(:,2),'.-')
coords0 = chain.coords;
coordsprev = coords0;
%%
for fc = opt.frames
        
    %% load in image and rescale
    img = obj.read(fc); 
    img = double(img(:,:,opt.col));
    img(:) = (img(:)-min(img(:)))/(max(img(:))-min(img(:)));

    %% adjust image brightness
    %img2 = imgaussfilt(img,1,'FilterSize',3);
    if (~isempty(opt.adjscl))
        img2 = imadjust(img,opt.adjscl,[0,1],opt.gamma);
    else
        img2 = img;
    end
    %imshow(img3,[],'InitialMagnification','fit');

    %% crop image
    if (~isempty(opt.croprect))
        img3 = imcrop(img2,opt.croprect);
    else
        img3 = img2;
    end
    
    %% get GVF field for image
    [X,Y] = meshgrid(1:size(img3,2),1:size(img3,1));
    [FextGVF,Fext] = getGVFfield(img3);
    
    if (opt.displaybeforeopt)
        imshow(img3,[],'InitialMagnification','fit');
        title(sprintf('Before Opt: Frame %d',fc))
        hold all
        
        if (chain.fixpos0)
            plot([chain.pos0(1); coords0(1:2:end); chain.posf(1)],[chain.pos0(2); coords0(2:2:end); chain.posf(2)],'.-')
            plot([chain.pos0(1); chain.coords(1:2:end); chain.posf(1)],[chain.pos0(2); chain.coords(2:2:end); chain.posf(2)],'.-')
        else
          plot( coordsprev(1:2:end), coordsprev(2:2:end),'.-')
            plot(chain.coords(1:2:end), chain.coords(2:2:end),'.-')
        end
        hold off
        drawnow
    end

    
    %% optimize filament
    optopt = opt;
    if (fc == opt.frames(1))
        optopt.displayevery = 100;
    end
    newcoords = optimizeFilament(chain,X,Y,FextGVF,img3,optopt)
    savecoords(:,fc) = newcoords;
    chain.coords = newcoords;
   
    %%
    if (opt.displayafteropt)
        imshow(img3,[],'InitialMagnification','fit');
        title(sprintf('Frame %d',fc))
        hold all
        if (chain.fixpos0)
            plot([chain.pos0(1); coords0(1:2:end); chain.posf(1)],[chain.pos0(2); coords0(2:2:end); chain.posf(2)],'.-')
            plot([chain.pos0(1); chain.coords(1:2:end); chain.posf(1)],[chain.pos0(2); chain.coords(2:2:end); chain.posf(2)],'.-')
        else
            plot(coordsprev(1:2:end), coordsprev(2:2:end),'.-')
             plot(chain.coords(1:2:end), chain.coords(2:2:end),'.-')
        end
        hold off
        drawnow
    end
    
     coordsprev = newcoords;

end
