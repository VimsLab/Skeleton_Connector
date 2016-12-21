function [ slices ] = readStack(filename, data_folder)
%Checks is MIP/filename exists in MIP_folder, if so, load it, if not, do a maximum
%intensity projection and return in
%Loads whole thing in memory 
%Also outputs slice with max intensity pixel

disp('Reading 3D Stack');

% if(exist([MIP_folder filename]))
%     mip_image = imread([MIP_folder filename]);
%     return;
% end

fname = [data_folder filename];
info = imfinfo(fname);
num_images = numel(info);
slice_one = imread(fname, 1, 'Info', info);
slices = zeros(size(slice_one,1), size(slice_one, 2), num_images);
slices(:,:,1) = slice_one;
for k = 2:num_images
    slices(:,:,k) = imread(fname, k, 'Info', info);
end
slices = uint8(slices);

end

