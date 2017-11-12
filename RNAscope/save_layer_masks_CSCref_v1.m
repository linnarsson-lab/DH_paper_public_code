tic 
close all
clear all


load ref_CSC
I = imread('mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/ref_CSC_colors_test.tif');
I = imresize(I,size(ref_CSC));
figure; imshow(I);
[n,x] = hist(double(I(:)),256);
[~,xi] = sort(n,'descend');
% top 8 colors are relvant
unival = sort(round(x(xi(1:8))));
bwall = false(size(ref_CSC));
for i=2:length(unival)
    bw = I==unival(i);
    if unival(i)==0
        bw = ~bw;
    end
    bw = imfill(bw,'holes');
    bw = bwareaopen(bw,3000);
    bw = imdilate(bw,strel('disk',11));
    layermask_bw{i} = bw;
    bwall = bwall +bw;
    layers_mask{i} = bw;
    figure; imshow(bw); hold on;
    B{i} = bwboundaries(bw);
    plot(B{i}{1}(:,2),B{i}{1}(:,1),'r');
    if i>2
        plot(B{i}{2}(:,2),B{i}{2}(:,1),'r');
        B{i} = B{i}(1:2);
    else
        B{i} = B{i}(1);
    end
end


bwall = imdilate(bwall,strel('disk',5));
layermask_bw{1} = bwall;
B{1} = bwboundaries(bwall);
figure; imshow(bwall); hold on;
plot(B{1}{1}(:,2),B{1}{1}(:,1),'r');


figure; imshow(ref_CSC); hold on;
for i=1:8
    for j=1:length(B{i})
        plot(B{i}{j}(:,2),B{i}{j}(:,1),'r');
    end
end
plot(set_ref_CSC(:,1),set_ref_CSC(:,2),'+y','markersize',20);
plot(set_ref_CSC(:,1),set_ref_CSC(:,2),'or','markersize',20);



save layers_mask_CSC_ref layermask_bw B

toc











