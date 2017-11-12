tic 
close all
clear all


load ref_LSC
I = imread('mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/ref_LSC_colors_test.tif');
I = imresize(I,size(ref_LSC));
figure; imshow(I);
[n,x] = hist(double(I(:)),256);
[~,xi] = sort(n,'descend');
% top 8 colors are relvant
unival = sort(round(x(xi(1:8))),'descend');
bwall = false(size(ref_LSC));
Ilayers_LSC = uint8(zeros(size(bwall)));
for i=1:length(unival)-1
    bw = I==unival(i);
    if unival(i)==0
        bw = ~bw;
    end
    bw = imfill(bw,'holes');
    bw = bwareaopen(bw,3000);
    bw = imdilate(bw,strel('disk',11));
    layermask_bw_LSC{i} = bw;
    bwall = bwall +bw;
    Ilayers_LSC(bw>0) = i;
%     layers_mask{i} = bw;
    figure; imshow(bw); hold on;
    B_LSC{i} = bwboundaries(bw);
    
    plot(B_LSC{i}{1}(:,2),B_LSC{i}{1}(:,1),'r');
    for j=1:length(B_LSC{i})
        B_LSC{i}{j}(:,2) = smoothn(B_LSC{i}{j}(:,2), 10);
        B_LSC{i}{j}(:,1) = smoothn(B_LSC{i}{j}(:,1), 10);
        plot(B_LSC{i}{j}(:,2),B_LSC{i}{j}(:,1),'r');

    end
end

i = length(unival);
bwall = imdilate(bwall,strel('disk',5));
layermask_bw_LSC{i} = bwall;
Ilayers_LSC(bwall==0) = i;
B_LSC{i} = bwboundaries(bwall);
B_LSC{i}{1}(:,2) = smoothn(B_LSC{i}{1}(:,2), 10);
B_LSC{i}{1}(:,1) = smoothn(B_LSC{i}{1}(:,1), 10);
figure; imshow(bwall); hold on;
plot(B_LSC{i}{1}(:,2),B_LSC{i}{1}(:,1),'r');


figure; imshow(ref_LSC); hold on;
for i=1:8
    for j=1:length(B_LSC{i})
        plot(B_LSC{i}{j}(:,2),B_LSC{i}{j}(:,1),'r');
    end
end
plot(set_ref_LSC(:,1),set_ref_LSC(:,2),'+y','markersize',20);
plot(set_ref_LSC(:,1),set_ref_LSC(:,2),'or','markersize',20);


save layers_mask_LSC_ref layermask_bw_LSC B_LSC Ilayers_LSC

figure('position',[10,10,1000,800],'color','w');
imagesc(Ilayers_LSC);
axis tight
axis equal
axis off


toc











