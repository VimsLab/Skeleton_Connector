function bw_conn = Connect_MST_ND(bw, gap_length, scale, endpointsOnly, window_sz, verbose)
% An implementation of the algorithm in [my paper citation] that connects
% nearby points in an N-dimensional skeleton or edge map. Useful if your
% segmentation algorithm contained gaps.
%
%Implemented in Matlab R2016b
%
%Input: 
%       bw - a logical N-dimensional array bw representing a skeleton or edge
%       points to be connected
%
%       gap_length (optional) - a single number representing 
%       the furthest gap the algorithm will connect. Default 50
%           
%       scale (optional) - a 1xN array, where N is the dimension of bw
%       containing the scale in each dimension. Useful for medical or
%       biological data where the scan resolution is different in the z
%       dimension. Default all 1s
%
%       endpointsOnly (optional) - a (nonlogical) flag. If the flag is set to 1,
%       only endpoints can be connected -- fast. If set to 0, the algorithm will
%       connect any points -- slower, especially on large skeletons.
%       If the flag is set to 2, connections between points must include at 
%       at least 1 endpoint -- also slower. Default 2
%
%       window_sz (optional) - a 1xN array, where N is the dimension of bw, 
%       containing a window size of how to break up the original bw array.
%       Useful if the original bw array is very large in memory.
%       Set all values <= 1 to run on the entire array. Default is all 30s
%
%       verbose (optional) - a logical flag (0 or 1) that determines if the
%       program prints to the screen its progress. Default 1
%
%Output: 
%       bw_conn - bw with added connections. 

%Check input to see if we have to set defaults.
if(nargin <= 1)
    gap_length = 50;
end
if(nargin <= 2)
    scale = ones(1, length(size(bw)));
end
if(nargin <= 3)
    endpointsOnly = 2;
end
if(nargin <= 4)
    window_sz = ones(1, length(size(bw)))*30;
end
if(nargin <= 5)
    verbose = 1;
end

if(all(window_sz <= 1))
    %If the window size is very small, simply run on the entire array
    if(verbose)
       disp('Running algorithm on entire array...'); 
    end
    bw_conn = Connect_MST_ND_Helper(bw,gap_length, scale, endpointsOnly);
else
    %if your array is very large, set window_sz to break it up 
    %into smaller, easier to process chunks that will fit in memory
    if(verbose)
       disp('Running algorithm on array chunks...'); 
    end
    outCell=mat2tiles(bw,window_sz);
    numCells = numel(outCell);
    for ii=1:numCells
        if(verbose)
            disp(['    Running algorithm on chunk: ' num2str(ii) '/' num2str(numCells)]);
        end
        if length(find(outCell{ii})) > 5
            outCell{ii} = Connect_MST_ND_Helper(outCell{ii},gap_length, scale, endpointsOnly);
        end
    end
    if(verbose)
        disp('Converting resulting cell back to array...');
    end
    bw_conn = cell2mat(outCell);
end

%Helper function that actually does the work for MST connections
function bw_conn = Connect_MST_ND_Helper(bw,gap_length, scale, endpointsOnly)

if(endpointsOnly == 1)
    %Endpoints connections
    [ep, l] = findEndpointsND(bw);
elseif(endpointsOnly == 0)
    %Connect to all points
    %Takes much more time (minutes vs seconds)
    ep = bw; l = find(ep);
elseif(endpointsOnly == 2)
    %Connections must include at least 1 endpoint.
    %Connect all points, but also identify true endpoints
    %Later, in the code, set weights between non-endpoints that are not connected
    %to 0
    ep = bw; l = find(ep);
    [ep2, l2] = findEndpointsND(bw);
    t = double(~ep2(l));
    %tM is square symmetrical matrix representing if connections are between
    %true endpoints or nonendpoints.
    tM = logical(t * t');
end



%Get list of indices of points we would like to connect
nd=length(size(bw));
subidx=cell(1,nd);
[subidx{:}] = ind2sub(size(bw), l);
idxs = [subidx{:}];
idxs_r_c_swapped = idxs;
idxs_r_c_swapped(:,[1 2]) = idxs_r_c_swapped(:,[2 1]);

%Create a graph that connects the points together
%Weights are euclidean distance
%The scale factor is used on each dimension
G = squareform(pdist(idxs_r_c_swapped.*repmat(scale, size(idxs_r_c_swapped, 1), 1)));

%For points that are already connected, make edge weight very small (eps)
CC = bwconncomp(bw);
for ii = 1:length(CC.PixelIdxList)
    pixels = CC.PixelIdxList{ii};
    [C, ia, ib] = intersect(l, pixels);
    G(ia,ia) = eps;
end
if(endpointsOnly == 2)
    %For non-endpoints that are not connected, break the edge connection
    G(tM & G~=eps)=0;
end

%Run the minimum spanning tree algorithm
[Tree, ~] = graphminspantree(sparse(G));
M = full(Tree);
bw_conn = bw;

%Valid connections are those that arent already connected and are less than
%gap length apart
[R,C] = find(M>2*eps & M<gap_length);
connections = [idxs(R,:), idxs(C,:)];

%Rasterize the points back onto the array
bw_conn = bw;
for ii = 1:length(R)
    %note connections is row-major -- r,c,dim3,dim4,dim5,...
    raster_pts = bresenhamND(connections(ii,1:nd), connections(ii,nd+1:end));
    lin_raster_pts =sub2indND(size(bw), raster_pts);
    bw_conn(lin_raster_pts) = 1;
end

