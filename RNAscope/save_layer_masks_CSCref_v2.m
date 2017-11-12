tic 
close all
clear all


load ref_CSC
I = imread('mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/ref_CSC_colors_test_v2.tif');
I = imresize(I,size(ref_CSC));
figure; imshow(I);
[n,x] = hist(double(I(:)),256);
[~,xi] = sort(n,'descend');
% top 8 colors are relvant
unival = sort(round(x(xi(1:8))),'descend');
bwall = false(size(ref_CSC));
Ilayers_CSC = uint8(zeros(size(bwall)));
for i=1:length(unival)-1
    bw = I==unival(i);
    if unival(i)==0
        bw = ~bw;
    end
    bw = imfill(bw,'holes');
    bw = bwareaopen(bw,3000);
    bw = imdilate(bw,strel('disk',11));
    layermask_bw_CSC{i} = bw;
    bwall = bwall +bw;
    Ilayers_CSC(bw>0) = i;
%     layers_mask{i} = bw;
    figure; imshow(bw); hold on;
    B_CSC{i} = bwboundaries(bw);
    
    plot(B_CSC{i}{1}(:,2),B_CSC{i}{1}(:,1),'r');
    for j=1:length(B_CSC{i})
        B_CSC{i}{j}(:,2) = smoothn(B_CSC{i}{j}(:,2), 10);
        B_CSC{i}{j}(:,1) = smoothn(B_CSC{i}{j}(:,1), 10);
        plot(B_CSC{i}{j}(:,2),B_CSC{i}{j}(:,1),'r');

    end
end

i = length(unival);
bwall = imdilate(bwall,strel('disk',5));
layermask_bw_CSC{i} = bwall;
Ilayers_CSC(bwall==0) = i;
B_CSC{i} = bwboundaries(bwall);
B_CSC{i}{1}(:,2) = smoothn(B_CSC{i}{1}(:,2), 10);
B_CSC{i}{1}(:,1) = smoothn(B_CSC{i}{1}(:,1), 10);
figure; imshow(bwall); hold on;
plot(B_CSC{i}{1}(:,2),B_CSC{i}{1}(:,1),'r');


figure; imshow(ref_CSC); hold on;
for i=1:8
    for j=1:length(B_CSC{i})
        plot(B_CSC{i}{j}(:,2),B_CSC{i}{j}(:,1),'r');
    end
end
plot(set_ref_CSC(:,1),set_ref_CSC(:,2),'+y','markersize',20);
plot(set_ref_CSC(:,1),set_ref_CSC(:,2),'or','markersize',20);


save layers_mask_CSC_ref layermask_bw_CSC B_CSC Ilayers_CSC

figure('position',[10,10,1000,800],'color','w');
imagesc(Ilayers_CSC);
axis tight
axis equal
axis off


toc
