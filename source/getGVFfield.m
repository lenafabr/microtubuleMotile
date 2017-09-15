function [FextGVF,Fext] = getGVFfield(img,options)
%% get GVF field for image
opt.Wline = 1;
opt.Wedge = 0;
opt.Wterm = 0;
opt.Sigma1 = 1;
opt.Sigma2 = 1;
opt.Mu = 0.2;
opt.Iterations = 100;
opt.getgvf = 1;

if (nargin > 1)
    opt = copyStruct(options,opt);
end
opt.SigmaGVF = 1;


% get gradient of energy function
Eext = ExternalForceImage2D(img,opt.Wline, opt.Wedge, opt.Wterm,opt.Sigma1);

Fx=ImageDerivatives2D(Eext,opt.Sigma2,'x');
Fy=ImageDerivatives2D(Eext,opt.Sigma2,'y');
Fext(:,:,1)=Fx*2*opt.Sigma2^2;
Fext(:,:,2)=Fy*2*opt.Sigma2^2;

% calculate the GVF field
if (opt.getgvf)
    FextGVF=GVFOptimizeImageForces2D(Fext, opt.Mu, opt.Iterations, opt.SigmaGVF);
else
    FextGVF = Fext;
end


end