function [compare] = pathEvaluate(intensity, position)
%UNTITLED6 Summary of this function goes here
%   calculate energy using bending energy function
%   determine good ratio between energy and brightness
%   output is a certain value that can be compared (higher is more
%   favorable) - (lower is less favorable)
%   Detailed explanation goes here

energy = BendingEnergy(position);

compare = 0.00081*intensity - energy; %play around with this output

end

