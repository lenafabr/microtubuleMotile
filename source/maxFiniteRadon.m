function [maxR,maxth,Rvals,thetavals] = maxFiniteRadon(img,pstart,L,ntheta)
% get maximal magnitude of finite length radon transform
% total length L, centered on point pstart
% ntheta: number of angles to sample; default=180;


% angles to sample (in degrees)
if (nargin>3)
    thetavals = linspace(0,180,ntheta+1);
    thetavals = thetavals(1:end-1);
else
    thetavals = 0:179;
end

% lengths to sample
lams = -L/2:1:L/2;

Rvals = zeros(size(thetavals));
for tc = 1:length(thetavals)
    theta = thetavals(tc)*pi/180;
    dvec = [cos(theta),sin(theta)];% direction
    linepts = [dvec(1)*lams'+pstart(1),dvec(2)*lams'+pstart(2)];
    Rvals(tc) = sum(interp2(img,linepts(:,1),linepts(:,2)));
end

[maxR,ind] = max(Rvals);
maxth = thetavals(ind);
%plot(thetavals,Rvals)