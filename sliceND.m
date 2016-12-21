function out = sliceND(A, ix, dim)

subses = repmat({':'}, [1 ndims(A)]);
subses{dim} = ix;
out = A(subses{:});
end
