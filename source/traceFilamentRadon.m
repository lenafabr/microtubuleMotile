function [filament, opt] = traceFilamentRadon(img,pstart,options)
% trace a filament using radon tranforms, from a given starting point
% all lengths are in pixes

% set default parameters
opt = struct();
opt.Lradon = 40; % length to use for finite radon transform
opt.Ltrace = 10; % separation between points for the tracing
% cutoff for radon transform magnitude to call filament end
% relative to starting point magnitude
opt.relRcutoff = 0.5; 
% cutoff in theta change to call filament end
% in degrees
opt.relthcutoff = 60;

% max points to trace in each direction
opt.maxpts = 100;

% display tracing during calculation
opt.dodisplay = 1;


%%
if(nargin>1); opt = copyStruct(options,opt); end
%%

if (opt.dodisplay)
    imshow(img,[],'InitialMagnification','fit');
end

% initial transform
[maxR0,maxth0,Rvals0,thvals0] = maxFiniteRadon(img, pstart,opt.Lradon);
theta = maxth0*pi/180;
dvec0 = [cos(theta),sin(theta)];% direction
    
if (opt.dodisplay)    
    lams = -opt.Lradon/2:1:opt.Lradon/2;
    linepts = [dvec0(1)*lams'+pstart(1),dvec0(2)*lams'+pstart(2)];
    hold all
    plot(pstart(1),pstart(2),'y.',linepts(:,1),linepts(:,2),'c')
    hold off
end
%
ptprev = pstart;
thprev = maxth0;
dvec = dvec0;

%% trace forward
maxRvals = [];
for tc = 1:opt.maxpts
    pt = ptprev + dvec*opt.Ltrace;
    
    [maxR,maxth] = maxFiniteRadon(img, pt,opt.Lradon);
    
    dth = maxth-thprev;
    if dth>90
        maxth = maxth-180;
    elseif dth<-90
        maxth = maxth+180;        
    end
    dth = maxth-thprev;    
    
    if (abs(dth)>opt.relthcutoff)
        % call filament end because too much reorientation
        break
    end
    
    maxRvals(tc) = maxR;
    
    theta = maxth*pi/180;
    dvec = [cos(theta),sin(theta)];% direction
    
    if (opt.dodisplay)
        lams = 0:1:opt.Ltrace;
        linepts = [dvec(1)*lams'+pt(1),dvec(2)*lams'+pt(2)];
        hold all
        plot(pt(1),pt(2),'y.',linepts(:,1),linepts(:,2),'g')
        hold off
        drawnow
    end
    
    ptprev = pt;
    thprev = maxth;
end

%% trace back
maxRvalsB = [];
ptprev=pstart;
thprev = maxth0;
dvec = dvec0;
for tc = 1:opt.maxpts
    pt = ptprev - dvec*opt.Ltrace;
    
    [maxR,maxth] = maxFiniteRadon(img, pt,opt.Lradon);
    
    dth = maxth-thprev;
    if dth>90
        maxth = maxth-180;
    elseif dth<-90
        maxth = maxth+180;        
    end
    dth = maxth-thprev;    
    
    if (abs(dth)>opt.relthcutoff)
        % call filament end because too much reorientation
        break
    end
    
    maxRvalsB(tc) = maxR;
    
    theta = maxth*pi/180;
    dvec = [cos(theta),sin(theta)];% direction
    
    if (opt.dodisplay)
        lams = -opt.Ltrace:1:0;
        linepts = [dvec(1)*lams'+pt(1),dvec(2)*lams'+pt(2)];
        hold all
        plot(pt(1),pt(2),'y.',linepts(:,1),linepts(:,2),'g')
        hold off
        drawnow
    end
    
    ptprev = pt;
    thprev = maxth;
end

end