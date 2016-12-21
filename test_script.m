clear; close all; clc;

%Read the 3D stack containing a skeleton
pth = './example_stack/';
img_name = 'example_skel.tif';
stack = readStack(img_name,pth);

%Set up parameters
gap_length = 100;
endpointsOnly = 2;
scale = [2.6 2.6 1.2]; %x y z
window_sz = round(scale*100);
Verbose = 1;

%Cal the MST algorithm on the skeleton
skeleton = logical(stack);
skeleton = bwareaopen(skeleton,5);
tic;
bw_conn = Connect_MST_ND(skeleton, gap_length, scale, endpointsOnly, window_sz, Verbose);
t = toc();
disp(['The MST algorithm took ' num2str(t) 's to process the ' num2str(size(skeleton)) ' size array']);


%Display results
if(length(size(skeleton)) == 3)
    figure;
    l = find(skeleton); [r c z] = ind2sub(size(skeleton), l); plot3(r*2.6, c*2.6, z*1.2, 'r.');
    hold on;
    l2 = find(bw_conn & ~skeleton); [r c z] = ind2sub(size(bw_conn), l2);plot3(r*2.6, c*2.6, z*1.2, 'b.');
    ax = gca; 
    ax.Clipping = 'off';
else
    figure;
    l = find(skeleton); [r c] = ind2sub(size(skeleton), l); plot(r*2.6, c*2.6, 'r.');
    hold on;
    l2 = find(bw_conn & ~skeleton); [r c] = ind2sub(size(bw_conn), l2);plot(r*2.6, c*2.6, 'b.');
    ax = gca; 
    ax.Clipping = 'off';
end
