%%
obj = VideoReader('../tubA_GFP_2mic_linear_cap15.mp4')

%% load in image
img = read(obj,5); 
%imshow(img,[],'InitialMagnification','fit');

%img = imread('frame1.png');
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

%% Initial snake points

ptinit = [221.91       234.03
       177.13       186.49
       183.76       159.95
       196.48       140.05];
  
%%
ptinit = getline

%%   
imshow(img3,[],'InitialMagnification','fit');
hold all
plot(ptinit(:,1),ptinit(:,2),'r.-')
hold off

% parameterize by arc length
param = arclenparam(ptinit');

%%
coords = interp1(1:size(ptinit,1),ptinit,linspace(1,size(ptinit,1),snake.npts),'spline')

%% set up chain object

chain.nseg = 20;
chain.nbead = chain.nseg-1;
chain.ncrd = 2*chain.nbead;

chain.len = param(end); % overall length of chain
chain.ls = chain.len/chain.nseg; % ground state segment length

% energetics
chain.lp = 1; % persistence length (nm)
chain.lstretch = 1; % stretch modulus (kT/nm)

% constraints
% fix end positions
chain.pos0 = ptinit(1,:)'; 
chain.posf = ptinit(end,:)';
% do not fix end tangents
chain.fixtan0 = 0; 
chain.fixtanf = 0;
chain.tan0 = [1,0];
chain.tanf = [1,0];

%% interpolate initial snake
% parameterize by arc length
param = arclenparam(ptinit');
% interpolate initial snake points
paramint = linspace(param(1),param(end),chain.nbead+2);
pts = interp1(param,ptinit,paramint,'linear')
chain.coords = reshape(pts(2:end-1,:)',2*chain.nbead,1)
hold all
plot(pts(:,1),pts(:,2),'.-')

coords0 = chain.coords
%% get energy and gradient for the interpolated chain
[energy,grad] = energyWLC2d(chain);

%% test gradient
[energy0,grad0] = energyWLC2d(chain);
tiny = 1e-6;
for c = 1:length(grad)
    chain.coords(c) = chain.coords(c) + tiny;
    [energy,grad] = energyWLC2d(chain);
    chain.coords(c) = chain.coords(c) - tiny;
    [c (energy-energy0)/tiny grad(c)]
end

%% gradient descent based on internal + external energies
stepsize = 0.2;
nstep = 10000;
%coords0 = chain.coords;
displayevery = 100;
extscl = 0.5;

% tolerance in rms step size per bead
% to define convergence
steptol = 1e-5;

for step = 1:nstep
    
    % gradient from internal energy
    [energy,grad] = energyWLC2d(chain);
    
    % GVF field
    FextX = interp2(X,Y,FextGVF(:,:,1),chain.coords(1:2:end),chain.coords(2:2:end));
    FextY = interp2(X,Y,FextGVF(:,:,2),chain.coords(1:2:end),chain.coords(2:2:end));
    %FextX = interp2(X,Y,Fext(:,:,1),chain.coords(1:2:end),chain.coords(2:2:end));
    %FextY = interp2(X,Y,Fext(:,:,2),chain.coords(1:2:end),chain.coords(2:2:end));
    
    gradext = reshape([FextY'; FextX'], 2*chain.nbead,1);
    
    % total internal + external gradient
    Ftot = -grad + gradext*extscl;
    
    rmsstep = norm(Ftot)/chain.nbead/2*stepsize;
    
    
    if (mod(step,displayevery)==0)
        disp([step rmsstep])
    end
    
    chain.coords = chain.coords + stepsize*Ftot;
    
    if (mod(step,displayevery) == 0)
    imshow(img3,[],'InitialMagnification','fit');
    hold all
    plot([chain.pos0(1); coords0(1:2:end); chain.posf(1)],[chain.pos0(2); coords0(2:2:end); chain.posf(2)],'.-')
    plot([chain.pos0(1); chain.coords(1:2:end); chain.posf(1)],[chain.pos0(2); chain.coords(2:2:end); chain.posf(2)],'.-')
    hold off
    drawnow
    end
        
    if (rmsstep<steptol)
        break
    end
end

%% test optimization as a function
options.steptol = 1e-4;
newcoords = optimizeFilament(chain,X,Y,FextGVF,img3,options)

%% view external gradient

% gradient from internal energy
    [energy,grad] = energyWLC2d(chain);
    
    % GVF field
    FextX = interp2(X,Y,FextGVF(:,:,1),chain.coords(1:2:end),chain.coords(2:2:end));
    FextY = interp2(X,Y,FextGVF(:,:,2),chain.coords(1:2:end),chain.coords(2:2:end));
    
    gradext = reshape([FextY'; FextX'], 2*chain.nbead,1);
    
    %% total internal + external gradient
    Ftot = -grad + gradext*extscl;
   
     imshow(img3,[],'InitialMagnification','fit');
    hold all
    plot(chain.coords(1:2:end),chain.coords(2:2:end),'.-')
    quiver(chain.coords(1:2:end),chain.coords(2:2:end),gradext(1:2:end),gradext(2:2:end))
    hold off
%% get GVF field for image
Options = struct();
Options.Wline = 1;
Options.Wedge = 0;
Options.Wterm = 0;
Options.Sigma1 = 1;
Options.Sigma2 = 1;
% get gradient of energy function
Eext = ExternalForceImage2D(img3,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1);

Fx=ImageDerivatives2D(Eext,Options.Sigma2,'x');
Fy=ImageDerivatives2D(Eext,Options.Sigma2,'y');
Fext(:,:,1)=Fx*2*Options.Sigma2^2;
Fext(:,:,2)=Fy*2*Options.Sigma2^2;

%% plot gradient field
[X,Y] = meshgrid(1:size(img3,2),1:size(img3,1));
%F1 = Fext(:,:,1); F2 = Fext(:,:,2);
imshow(Eext,[])
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
quiver(X,Y,Fext(:,:,2),Fext(:,:,1),4)
hold off
