function [valid_AngPeaks] = RadonTransform(img,pstart,L, prev_angle, opt)
% get maximal magnitude of finite length radon transform
% total length L, centered on point pstart
% ntheta: number of angles to sample; default=180;

thEnergy = opt.thEnergy;
thRelEnergy = opt.thRelEnergy;
thAngle = opt.thAngle;% was 70
thAngleDiff = opt.thAngleDiff;
nPeaks = opt.nPeaks;
%change thetavals for angle cut off?
thetavals = 0:359;




% lengths to sample
lams = 0:1:L;

Rvals = zeros(size(thetavals));
for tc = 1:length(thetavals)
    theta = thetavals(tc)*pi/180;
    dvec = [cos(theta),sin(theta)];% direction
    linepts = [dvec(1)*lams'+pstart(1),dvec(2)*lams'+pstart(2)];
    Rvals(tc) = sum(interp2(img,linepts(:,1),linepts(:,2)));
end





%find peaks and threshold peaks to 
peaks = findpeaks(Rvals);
plength = size(peaks',1);

tangles = zeros;

%example  
for i = 1:plength
    index = find(Rvals==peaks(i)); %****resolve issue here****
    tvidx = thetavals(index(1)); 
    tangles(i) = tvidx;
end
    
peaks = [tangles; peaks];

[valid_AngPeaks] =  peakValidation(thRelEnergy, thEnergy, thAngle, thAngleDiff, peaks, nPeaks,prev_angle );
% See graphically
%disp(tangles)
%plot(tangles, peaks, 'or')
