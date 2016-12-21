function [ ep, lin_ind ] = findEndpointsND( bw, edgeEndpoints )
%Function to quickly find all endpoints of a logical N-dimensional array
%containing a skeleton
%
%Input:
%       bw - a logical N-dimensional array containing a skeleton (thin,
%       only is long in 1 dimension
%
%       edgeEndpoints (optional) - a logical flag that determines if lines
%       ending on the edge of the array count as endpoints. If 1, these
%       will count, if 0, the array is padded with 1s so the edges of the array
%       will not be endpoints. Default is 0
%
%Output: 
%       ep - an array the same size and type as bw, with a 1 if there is an
%       endpoint at the location.
%
%       lin_ind - linear indices of the endpoint locations corresponding to
%       bw

if(nargin <= 1)
    edgeEndpoints = 0;
end

%Create a filter that is 3x3x3x3x....
%With a 1 at all locations, and a 10 at the center
n = ndims(bw);
shape3s = ones(1,n)*3;
filt = ones(shape3s); 
center_lin = round(length(filt(:)) / 2);
nms = numel(filt)+1; %Make the center 1 larger than the number of elements
filt(center_lin) = nms;

%Convolve the filter with the bw array
%If both the center and 1 side point is hit (i.e nms+1), then
% that is an endpoint
%Pad bw first, so that a line at the edge of the window doesnt get counted
%as an endpoint
if(~edgeEndpoints)
    bw = padarray(bw, ones(1,n), 1);
    ep = ismember(convn(double(bw),filt,'same'),nms+1);
    idx = arrayfun(@(sz)2:sz-1,size(ep),'UniformOutput',false);
    ep = ep(idx{:});
else
     ep = ismember(convn(double(bw),filt,'same'),nms+1);
end

lin_ind = find(ep);

end

