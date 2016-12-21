function out_pts = bresenhamND(pt1,pt2)
%Vectorized function to rasterize a line between 2 end points in N-dimensions.
%Input: 
%       pt1 - Inital point 1 X N vector.
%       pt2 - End point 1 X N vector.
%Output: pts - positions to raterize line M X N matrix (M number of N-dimensional points). 
%           x1 y1 z1 w1 ...
%           x2 y2 z2 w2 ...
%           ... ...  ... ...

direction = pt2-pt1; %get direction.
num_steps = max(abs(direction),[],2); %find number of steps to take.
step_size = direction/num_steps; %scale the direction for step size.

%starting from "pt1" step along the "direction" for "num_steps" with "step_size" intervals. 
out_pts = round(bsxfun(@plus,pt1,(0:num_steps)'*step_size)); %pixel positions to rasterize line.