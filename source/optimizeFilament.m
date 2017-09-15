function newcoords = optimizeFilament(chain,X,Y,FextGVF,img,options)
% optimize a filament position in a single frame
% uses gradient vector flow

opt.stepsize = 0.2;
opt.nstep = 10000;
%coords0 = chain.coords;
opt.displayevery = 100;
opt.extscl = 0.5;

% tolerance in rms step size per bead
% to define convergence
opt.steptol = 1e-4;

if (nargin>3)
     opt = copyStruct(options, opt);
end

coords0 = chain.coords;

for step = 1:opt.nstep
    
    % gradient from internal energy
    [energy,grad] = energyWLC2d(chain);
    
    % GVF field
    FextX = interp2(X,Y,FextGVF(:,:,1),chain.coords(1:2:end),chain.coords(2:2:end));
    FextY = interp2(X,Y,FextGVF(:,:,2),chain.coords(1:2:end),chain.coords(2:2:end));
    %FextX = interp2(X,Y,Fext(:,:,1),chain.coords(1:2:end),chain.coords(2:2:end));
    %FextY = interp2(X,Y,Fext(:,:,2),chain.coords(1:2:end),chain.coords(2:2:end));
    
    gradext = reshape([FextY'; FextX'], 2*chain.nbead,1);
    
    % total internal + external gradient
    Ftot = -grad + gradext*opt.extscl;
    
    if (mod(step,opt.displayevery)==0)
        disp([step rmsstep])
    end
    
    % edge beads cannot move along chain contour (no reptation)
    if (~chain.fixpos0)
        v1 = chain.coords(1:2)-chain.coords(3:4);
        v1 = v1/norm(v1);
        Ftot(1:2) = Ftot(1:2) - (v1'*Ftot(1:2))*v1;
    end
    if (~chain.fixposf)
        v1 = chain.coords(end-1:end)-chain.coords(end-3:end-2);
        v1 = v1/norm(v1);
        Ftot(end-1:end) = Ftot(end-1:end) - (v1'*Ftot(end-1:end))*v1;
    end
    
    rmsstep = norm(Ftot)/chain.nbead/2*opt.stepsize;
    
    chain.coords = chain.coords + opt.stepsize*Ftot;
    
    if (mod(step,opt.displayevery) == 0)
    imshow(img,[],'InitialMagnification','fit');
    hold all
    if (chain.fixpos0)
        plot([chain.pos0(1); coords0(1:2:end); chain.posf(1)],[chain.pos0(2); coords0(2:2:end); chain.posf(2)],'.-')
        plot([chain.pos0(1); chain.coords(1:2:end); chain.posf(1)],[chain.pos0(2); chain.coords(2:2:end); chain.posf(2)],'.-')
    else
        plot( coords0(1:2:end), coords0(2:2:end),'.-')
        plot(chain.coords(1:2:end), chain.coords(2:2:end),'.-')
    end
    hold off
    drawnow
    end
        
    if (rmsstep<opt.steptol)
        break
    end
end

newcoords = chain.coords;
end