function lin_idx = sub2indND(sizeM, J )%
%sub2indND converts an array of indices J for an N-dimensional array M to
%linear indices
%
%Input: 
%       sizeM - 1xN array containing the sizes of each dimension
%       J - PxN matrix containing P points with N dimensions
%           r1 c1 z1 w1 ...
%           r2 c2 z2 w2 ...
%           ... ...  ....
%
%Output:
%   lin_idx - Px1 array of linear indices
%
%
% Implements the following equation in a vectorized manner:
%       lin = i1 + (i2-1)*size(A,1) + (i3-1)*size(A,1)*size(A,2) + ... 

J(:, 2:end) = J(:,2:end)-1;
lin_idx = sum(repmat(cumprod([1 sizeM(1:end-1)]), [size(J,1), 1]).*J, 2);

end