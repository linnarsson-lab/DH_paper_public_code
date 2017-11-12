function [] = RNAscope_quant_dorsalhorn_v2(path_folder,filename,chan1name,chan2name,chan3name,chan4name,manual_th,fig_flag)

% path_folder = '/data2/C1_stuff/Dorsal_horn_MH/RNAscope/cFos_in_SNI/Cfos_488_Gad2_555_Slc17a6_643_M6_Sham_VIII2_LSC3/';
% filename = 'MergedStitched_Cfos_488_Gad2_555_Slc17a6_643_M6_Sham_VIII2_LSC3_cropped.tif';
% chan1name = 'gad2';
% chan2name = 'slc17a6';
% chan3name = 'dapi';
% chan4name = 'fos';
% manual_th = 2;

nbit = 16;
dapi = imread([path_folder,filename],'index',3);
dapi = imrotate(dapi,-90);
% dapi = dapi(4630:end,1900:end);
low_in = 110/(2^nbit-1); %double(prctile(dapi(:),60))/(2^nbit-1);
high_in = 450/(2^nbit-1); %double(prctile(dapi(:),99.9))/(2^nbit-1);
dapi = imadjust(dapi,[low_in,high_in],[0,1]);

geneC1 = imread([path_folder,filename],'index',1);
geneC1 = imrotate(geneC1,-90);
% geneC1 = geneC1(4630:end,1900:end);
% low_in = 100/(2^nbit-1); %double(prctile(geneC1(:),90))/(2^nbit-1);
% high_in = 450/(2^nbit-1); %double(prctile(geneC1(:),100))/(2^nbit-1);
low_in = double(prctile(geneC1(:),0.1))/(2^nbit-1);
high_in = double(prctile(geneC1(:),100))/(2^nbit-1);
geneC1 = imadjust(geneC1,[low_in,high_in],[0,1]);
T = adaptthresh(geneC1, 0.05);
geneC1 = geneC1-uint16(T*2^16);
geneC1 = imadjust(geneC1,[0,double(prctile(geneC1(:),99.99))/(2^nbit-1)],[0,1]);
% T = uint16(adaptthresh(geneC1, 0.05)*2^16);
% geneC1 = geneC1 - imadjust(T,double([min(T(:)),max(T(:))])./(2^nbit-1),[0,1]);
% geneC1 = medfilt2(geneC1);
% geneC1 = imadjust(geneC1,[0,double(prctile(geneC1(:),99.99))/(2^nbit-1)],[0,1]);

geneC2 = imread([path_folder,filename],'index',2);
geneC2 = imrotate(geneC2,-90);
% geneC2 = geneC2(4630:end,1900:end);
% low_in = 100/(2^nbit-1); %double(prctile(geneC2(:),90))/(2^nbit-1);
% high_in = 250/(2^nbit-1); %double(prctile(geneC2(:),100))/(2^nbit-1);
low_in = double(prctile(geneC2(:),0.1))/(2^nbit-1);
high_in = double(prctile(geneC2(:),100))/(2^nbit-1);
geneC2 = imadjust(geneC2,[low_in,high_in],[0,1]);
T = adaptthresh(geneC2, 0.05);
geneC2 = geneC2-uint16(T*2^16);
geneC2 = imadjust(geneC2,[0,double(prctile(geneC2(:),99.99))/(2^nbit-1)],[0,1]);
% geneC2 = geneC2 - imadjust(uint16(adaptthresh(geneC2, 0.05)*2^16));
% geneC2 = medfilt2(geneC2);
% geneC2 = imadjust(geneC2,[0,double(prctile(geneC2(:),99.99))/(2^nbit-1)],[0,1]);


geneC4 = imread([path_folder,filename],'index',4);
geneC4 = imrotate(geneC4,-90);
% geneC4 = geneC4(4630:end,1900:end);
% low_in = 100/(2^nbit-1); %double(prctile(geneC4(:),98))/(2^nbit-1);
% high_in = 450/(2^nbit-1); %double(prctile(geneC4(:),100))/(2^nbit-1);
low_in = double(prctile(geneC4(:),0.1))/(2^nbit-1);
high_in = double(prctile(geneC4(:),100))/(2^nbit-1);
geneC4 = imadjust(geneC4,[low_in,high_in],[0,1]);
T = adaptthresh(geneC4, 0.05);
geneC4 = geneC4-uint16(T*2^16);
geneC4 = imadjust(geneC4,[0,double(prctile(geneC4(:),99.99))/(2^nbit-1)],[0,1]);
% geneC4 = geneC4 - imadjust(uint16(adaptthresh(geneC4, 0.05)*2^16));
% geneC4 = medfilt2(geneC4);
% geneC4 = imadjust(geneC4,[0,double(prctile(geneC4(:),99.99))/(2^nbit-1)],[0,1]);
% % % % % % % % % % % % % % % % % % % % % % % 

I_rgb(:,:,1) = geneC1+dapi/3;
I_rgb(:,:,2) = geneC4+dapi/3;
I_rgb(:,:,3) = geneC2+dapi/3;

bwgeneC1 = imextendedmax(geneC1,20000, 4);
bwgeneC1 = imfill(bwgeneC1,'holes');
bwgeneC1 = bwgeneC1 - bwareaopen(bwgeneC1,50);

bwgeneC2 = imextendedmax(geneC2,10000, 4);
bwgeneC2 = imfill(bwgeneC2,'holes');
bwgeneC2 = bwgeneC2 - bwareaopen(bwgeneC2,50);

bwgeneC4 = imextendedmax(geneC4,20000, 4);
bwgeneC4 = imfill(bwgeneC4,'holes');
bwgeneC4 = bwgeneC4 - bwareaopen(bwgeneC4,50);


rgb_bw(:,:,1) = bwgeneC1;
rgb_bw(:,:,2) = bwgeneC4;
rgb_bw(:,:,3) = bwgeneC2;
rgb_bw = uint8(rgb_bw)*255;
figure('visible',fig_flag); imshow(rgb_bw); hold on;

rgb_bw_w(:,:,1) = ~bwgeneC4 & ~bwgeneC2;%~bwmaf;
rgb_bw_w(:,:,2) = ~bwgeneC1 & ~bwgeneC2;%~bwgeneC1;
rgb_bw_w(:,:,3) = ~bwgeneC1 & ~bwgeneC4;%~bwgeneC2;
rgb_bw_w = uint8(rgb_bw_w)*255;
figure('visible',fig_flag); imshow(rgb_bw_w); hold on;

geneC1_positions = regionprops(bwgeneC1>0,'Centroid','Area');
geneC1_positions = reshape([geneC1_positions.Centroid]',2,length(geneC1_positions))';
geneC4_positions = regionprops(bwgeneC4>0,'Centroid','Area');
geneC4_positions = reshape([geneC4_positions.Centroid]',2,length(geneC4_positions))';
geneC2_positions = regionprops(bwgeneC2>0,'Centroid','Area');
geneC2_positions = reshape([geneC2_positions.Centroid]',2,length(geneC2_positions))';
% % % % % % % % % % % % % % % % % % % % % % % % dapi segmentation
if isempty(dir([path_folder,filename(1:end-4),'_dapi_segmentation.mat']))
    bw = im2bw(dapi,double(prctile(dapi(:),70))/(2^nbit-1));%im2bw(dapi,graythresh(dapi));
    bw1 = bwareaopen(bw,10);
    bw1 = imfill(bw1,'holes');
    bw1 = bwareaopen(bw1,100);
    
    [bw2] = segment_dapi_by_gradient(dapi,bw1,100,8000,0);
    
    stats_dapi = regionprops(bw2,'area','boundingbox','PixelList','PixelIdxList');
    dapi_a = [stats_dapi.Area];
    bb_bw1 = [reshape([stats_dapi.BoundingBox]',4,length(stats_dapi))]';
    bw3 = bw2;
    for k=1:length(dapi_a)
        k
        if dapi_a(k)>8000
            box_ind_r = round([bb_bw1(k,2):(bb_bw1(k,2)+bb_bw1(k,4)-1)]');
            box_ind_c = round([bb_bw1(k,1):(bb_bw1(k,1)+bb_bw1(k,3)-1)]');
            threshold = prctile(dapi(stats_dapi(k).PixelIdxList),50);
            tmp = dapi(box_ind_r,box_ind_c);
            linInd = sub2ind(size(tmp),round(stats_dapi(k).PixelList(:,2)-bb_bw1(k,2)),round(stats_dapi(k).PixelList(:,1)-bb_bw1(k,1)));
            tmp(setdiff([1:length(tmp(:))],linInd)) = 0;
            BW4 = im2bw(tmp,double(threshold)/(2^nbit-1));
            BW4 = imfill(BW4,'holes');
            Dbw1 = -bwdist(~BW4); %Dbw1(~bw1) = -Inf;
            mask = imextendedmin(Dbw1,2);
            D2 = imimposemin(Dbw1,mask);
            L = watershed(D2);
            BW4(L == 0) = 0;
            bw3(stats_dapi(k).PixelIdxList) = BW4(linInd);
        end
    end
    
    
    save([path_folder,filename(1:end-4),'_dapi_segmentation'],'bw3')
else
    load([path_folder,filename(1:end-4),'_dapi_segmentation'])
end
% % % % % % % % % % % % % % % % % % % % % % % % % end dapi segmentation


stats = regionprops(bw3,'area');
dapi_a = [stats.Area];
B = bwboundaries(bw3);
figure('visible',fig_flag);
%     imshow(true(size(dapi))); hold on;
imshow(I_rgb); hold on;
for k = 1:length(B)
    boundary = B{k};
    if dapi_a(k)<100
        plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 0.5)
    elseif dapi_a(k)<15000
        B{k} = [smoothn(B{k}(:,1),50),smoothn(B{k}(:,2),50)];
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', 0.5)
    else
        plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 0.5)
    end
end
L = bwlabel(bw3);
stats_bw3 = regionprops(bw3,'Area','Image','BoundingBox','centroid');
bb_bw3 = [reshape([stats_bw3.BoundingBox]',4,length(stats_bw3))]';
dapi_cent = reshape([stats_bw3.Centroid]',2,length(stats_bw3))';
xmedline = length(dapi(1,:))/2;

geneC1_cell = zeros(1,max(L(:)));
geneC4_cell = zeros(1,max(L(:)));
geneC2_cell = zeros(1,max(L(:)));
geneC1_dots_cell = zeros(1,max(L(:)));
geneC4_dots_cell = zeros(1,max(L(:)));
geneC2_dots_cell = zeros(1,max(L(:)));
%     figure; imagesc(L);
for i=1:max(L(:));
    i
    box_ind_r = round([bb_bw3(i,2):(bb_bw3(i,2)+bb_bw3(i,4)-1)]');
    box_ind_c = round([bb_bw3(i,1):(bb_bw3(i,1)+bb_bw3(i,3)-1)]');
    area_cell(i) = stats_bw3(i).Area;
    if area_cell(i)<15000 & area_cell(i)>100
        tmp = geneC1(box_ind_r,box_ind_c);        
        geneC1_cell(i) = sum(tmp(stats_bw3(i).Image));%/area_cell(i);
        tmp = geneC4(box_ind_r,box_ind_c); 
        geneC4_cell(i) = sum(tmp(stats_bw3(i).Image));%/area_cell(i);
        tmp = geneC2(box_ind_r,box_ind_c); 
        geneC2_cell(i) = sum(tmp(stats_bw3(i).Image));%/area_cell(i);        
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%         tmp = regionprops(bwgeneC1(box_ind_r,box_ind_c)>0,'centroid');
        tmp = regionprops(bwgeneC1(box_ind_r,box_ind_c) & stats_bw3(i).Image,'centroid');
        geneC1_dots_cell(i) = length(tmp);%tmp.NumObjects;
%         if geneC1_dots_cell(i)>0
%             tmp = reshape([tmp.Centroid]',2,length(tmp))';
%             plot(tmp(:,1)+bb_bw3(i,1),tmp(:,2)+bb_bw3(i,2),'o','markeredgecolor','w','markerfacecolor','g');
%         end
%         tmp = regionprops(bwgeneC4(box_ind_r,box_ind_c)>0,'centroid');
        tmp = regionprops(bwgeneC4(box_ind_r,box_ind_c) & stats_bw3(i).Image,'centroid');
        geneC4_dots_cell(i) = length(tmp);%tmp.NumObjects;
%         if geneC4_dots_cell(i)>0
%             tmp = reshape([tmp.Centroid]',2,length(tmp))';
%             plot(tmp(:,1)+bb_bw3(i,1),tmp(:,2)+bb_bw3(i,2),'p','markeredgecolor','w','markerfacecolor','r');
%         end
%         tmp = regionprops(bwgeneC2(box_ind_r,box_ind_c)>0,'centroid');
        tmp = regionprops(bwgeneC2(box_ind_r,box_ind_c) & stats_bw3(i).Image,'centroid');
        geneC2_dots_cell(i) = length(tmp);%tmp.NumObjects;
%         if geneC2_dots_cell(i)>0
%             tmp = reshape([tmp.Centroid]',2,length(tmp))';
%             plot(tmp(:,1)+bb_bw3(i,1),tmp(:,2)+bb_bw3(i,2),'d','markeredgecolor','w','markerfacecolor','y');
%         end
    end
end


geneC1_all = [geneC1_cell];
geneC4_all = [geneC4_cell];
geneC2_all = [geneC2_cell];
geneC1_dots_all = [geneC1_dots_cell];
geneC4_dots_all = [geneC4_dots_cell];
geneC2_dots_all = [geneC2_dots_cell];


% pos_th = 0;

figure('visible',fig_flag);
set(gcf,'position',[100,100,1350,450],'color','w')
subplot(1,3,1);
tmp1 = geneC2_dots_all; 
tmp2 = geneC1_dots_all; 
tmp1_s = sort(tmp1(tmp1>0),'descend');
tmp2_s = sort(tmp2(tmp2>0),'descend');
plot(tmp1,tmp2,'o','markerfacecolor','r','markersize',3); hold on;
% [~,imin] = min( tmp1_s.^2 + ([1:length(tmp1_s)]).^2); 
thtmp1 = manual_th;%prctile(tmp1_s,pos_th);% tmp1_s(imin);
[~,imin] = min( tmp2_s.^2 + ([1:length(tmp2_s)]).^2); 
thtmp2 = manual_th;%prctile(tmp2_s,pos_th);%tmp2_s(imin);
plot(tmp1(tmp2>thtmp2 & tmp1>thtmp1),tmp2(tmp2>thtmp2 & tmp1>thtmp1),'sg','markerfacecolor','g');hold on;
plot(thtmp1*[1,1],[0,max(tmp2)],'k')
plot([0,max(tmp1)],thtmp2*[1,1],'k')
axis tight
xlabel([chan2name,' (dots/cell)']);ylabel([chan1name,' (dots/cell)']);
title([chan2name,' (pos)',num2str(sum(tmp1>thtmp1)),', ',' ',chan1name,'(pos)',num2str(sum(tmp2>thtmp2))...
    ,', both=',num2str(sum(tmp2>thtmp2 & tmp1>thtmp1))]);
ind_double_pos_geneC2_geneC1 = find(tmp2>thtmp2 & tmp1>thtmp1);
subplot(1,3,2);
tmp1 = geneC2_dots_all;
tmp2 = geneC4_dots_all;
tmp1_s = sort(tmp1(tmp1>0),'descend');
tmp2_s = sort(tmp2(tmp2>0),'descend');
plot(tmp1,tmp2,'o','markerfacecolor','r','markersize',3); hold on;
thtmp1 =manual_th;%prctile(tmp1_s,pos_th);%tmp1_s(imin);
thtmp2 = manual_th;%prctile(tmp2_s,pos_th);%tmp2_s(imin);
plot(tmp1(tmp2>thtmp2 & tmp1>thtmp1),tmp2(tmp2>thtmp2 & tmp1>thtmp1),'sg','markerfacecolor','g');hold on;
plot(thtmp1*[1,1],[0,max(tmp2)],'k')
plot([0,max(tmp1)],thtmp2*[1,1],'k')
axis tight
xlabel([chan2name,' (dots/cell)']);ylabel([chan4name,' (dots/cell)']);
title([chan2name,' (pos)',num2str(sum(tmp1>thtmp1)),', ',' ',chan4name,'(pos)'....
    ,num2str(sum(tmp2>thtmp2)),', both=',num2str(sum(tmp2>thtmp2 & tmp1>thtmp1))]);
ind_double_pos_geneC2_geneC4 = find(tmp2>thtmp2 & tmp1>thtmp1);
subplot(1,3,3);
tmp1 = geneC1_dots_all;
tmp2 = geneC4_dots_all;
tmp1_s = sort(tmp1(tmp1>0),'descend');
tmp2_s = sort(tmp2(tmp2>0),'descend');
plot(tmp1,tmp2,'o','markerfacecolor','r','markersize',3); hold on;
thtmp1 = manual_th;% prctile(tmp1_s,pos_th);%tmp1_s(imin);
thtmp2 = manual_th;%prctile(tmp2_s,pos_th);% tmp2_s(imin);
plot(tmp1(tmp2>thtmp2 & tmp1>thtmp1),tmp2(tmp2>thtmp2 & tmp1>thtmp1),'sg','markerfacecolor','g');hold on;
plot(thtmp1*[1,1],[0,max(tmp2)],'k')
plot([0,max(tmp1)],thtmp2*[1,1],'k')
axis tight
xlabel([chan1name,' (dots/cell)']);ylabel([chan4name,' (dots/cell)']);
title([chan1name,' (pos)',num2str(sum(tmp1>thtmp1)),', ',' ',chan4name,'(pos)'...
    ,num2str(sum(tmp2>thtmp2)),', both=',num2str(sum(tmp2>thtmp2 & tmp1>thtmp1))]);
ind_double_pos_geneC1_geneC4 = find(tmp2>thtmp2 & tmp1>thtmp1);

eval(['export_fig ',path_folder,filename(1:end-4),'scatter_quantification_dots_',date,'.pdf']);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
pos_cells_mat = [geneC1_dots_all'>manual_th, geneC2_dots_all'>manual_th, geneC4_dots_all'>manual_th];
[pos_cells_mat_sorted,XI] = sortrows(pos_cells_mat);
[uni_comb,ia] = unique(pos_cells_mat_sorted,'rows','first');
[~,ib] = unique(pos_cells_mat_sorted,'rows','last');

legend_str = repmat({''},length(uni_comb(:,1))-1,1);
for i=1:length(legend_str)
    if uni_comb(i+1,1)==1
        legend_str{i} = [legend_str{i},'&',chan1name];
    end
    if uni_comb(i+1,2)==1
        legend_str{i} = [legend_str{i},'&',chan2name];
    end
    if uni_comb(i+1,3)==1
        legend_str{i} = [legend_str{i},'&',chan4name];
    end
    legend_str{i} = legend_str{i}(2:end);
end

figure('visible',fig_flag);
imshow(I_rgb); hold on;
colorvec = distinguishable_colors(length(legend_str));
colorvec = flipud(colorvec);
for i=2:length(legend_str)+1
    indtmp = XI(ia(i):ib(i));
    for k = 1
        boundary = B{indtmp(k)};
        plot(boundary(:,2)+1, boundary(:,1)+1, 'color',colorvec(i-1,:), 'LineWidth', 0.5)
    end
end
for i=2:length(legend_str)+1
    indtmp = XI(ia(i):ib(i));
    for k = 1:length(indtmp)
        boundary = B{indtmp(k)};
        plot(boundary(:,2)+1, boundary(:,1)+1, 'color',colorvec(i-1,:), 'LineWidth', 0.5)
    end
end
legend(legend_str);
eval(['export_fig ',path_folder,filename(1:end-4),'_types_overlayed_RNAscope_',date,'.pdf']);


figure('color','w','position',[1,1,2000,1200],'visible',fig_flag);hold on;
colorvec = distinguishable_colors(length(legend_str));
colorvec = flipud(colorvec);
for i=2:length(legend_str)+1
    indtmp = XI(ia(i):ib(i));
    for k = 1
        boundary = B{indtmp(k)};
        plot(boundary(:,2)+1, boundary(:,1)+1, 'color',colorvec(i-1,:), 'LineWidth', 0.5)
    end
end
plot(geneC1_positions(:,1), geneC1_positions(:,2) ,'.r','markersize',0.1);hold on;
plot(geneC2_positions(:,1), geneC2_positions(:,2) ,'.b','markersize',0.1);hold on;
plot(geneC4_positions(:,1), geneC4_positions(:,2) ,'.g','markersize',0.1);hold on;
for i=2:length(legend_str)+1
    indtmp = XI(ia(i):ib(i));
    for k = 1:length(indtmp)
        boundary = B{indtmp(k)};
        plot(boundary(:,2)+1, boundary(:,1)+1, 'color',colorvec(i-1,:), 'LineWidth', 0.5)
    end
end

set(gca,'ydir','reverse','xlim',[0,length(dapi(1,:))],'ylim',[0,length(dapi(:,1))])
axis equal
axis off
legend(legend_str); 
eval(['export_fig ',path_folder,filename(1:end-4),'_types_whitebackground_',date,'.pdf']);



figure('color','w','position',[1,1,2000,1200],'visible',fig_flag);hold on;
plot(geneC1_positions(:,1), geneC1_positions(:,2) ,'.r','markersize',3);hold on;
plot(geneC2_positions(:,1), geneC2_positions(:,2) ,'.b','markersize',3);hold on;
plot(geneC4_positions(:,1), geneC4_positions(:,2) ,'.g','markersize',3);hold on;
for k = 1:length(B)
    boundary = B{k};
    if dapi_a(k)<100
%         plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 0.5)
    elseif dapi_a(k)<15000
        plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', 0.5)
    else
%         plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 0.5)
    end
end
legend(chan1name,chan2name,chan4name); 
set(gca,'ydir','reverse','xlim',[0,length(dapi(1,:))],'ylim',[0,length(dapi(:,1))])
axis equal
axis off
eval(['export_fig ',path_folder,filename(1:end-4),'_dapi_and_Dots_',date,'.pdf']);



figure('visible',fig_flag,'color','w','position',[1,1,2000,1200]);hold on;
colorvec = distinguishable_colors(length(legend_str));
colorvec = flipud(colorvec);
for i=2:length(legend_str)+1
    indtmp = XI(ia(i):ib(i));
    for k = 1
        boundary = B{indtmp(k)};
        plot(boundary(:,2)+1, boundary(:,1)+1, 'color',colorvec(i-1,:), 'LineWidth', 0.5)
    end
end
for i=2:length(legend_str)+1
    indtmp = XI(ia(i):ib(i));
    for k = 1:length(indtmp)
        boundary = B{indtmp(k)};
        plot(boundary(:,2)+1, boundary(:,1)+1, 'color',colorvec(i-1,:), 'LineWidth', 0.5)
    end
end
set(gca,'ydir','reverse','xlim',[0,length(dapi(1,:))],'ylim',[0,length(dapi(:,1))])
axis equal
axis off
legend(legend_str); 
eval(['export_fig ',path_folder,filename(1:end-4),'_types_noDots_whitebackground_',date,'.pdf']);


% colorvec = distinguishable_colors(length(legend_str));
% colorvec = flipud(colorvec);
% for i=2:length(legend_str)+1
%     figure('visible',fig_flag,'color','w','position',[1,1,2000,1200]);hold on;
%     indtmp = XI(ia(i):ib(i));
%     for k = 1:length(indtmp)
%         boundary = B{indtmp(k)};
%         plot(boundary(:,2)+1, boundary(:,1)+1, 'color',colorvec(i-1,:), 'LineWidth', 0.5)
%     end
%     set(gca,'ydir','reverse','xlim',[0,length(dapi(1,:))],'ylim',[0,length(dapi(:,1))])
%     axis equal
%     axis off
%     title(legend_str{i-1})
%     eval(['export_fig ',path_folder,filename(1:end-4),legend_str{i-1},'_whitebackground_',date,'.pdf']);
% end



sum_combination = zeros(1,length(legend_str));
for i=2:length(legend_str)+1
    indtmp = XI(ia(i):ib(i));
    sum_combination(i-1) = ib(i)-ia(i)+1;
end
table1 = [[{'dapi'};{sum(area_cell>100)}],[legend_str';m2c(sum_combination)]];
saveCellFile(table1,[path_folder,filename(1:end-4),'_counts_',date,'.txt']);

pos_cells_mat = [geneC1_dots_all(dapi_cent(:,1)<xmedline)'>manual_th, ....
     geneC2_dots_all(dapi_cent(:,1)<xmedline)'>manual_th, geneC4_dots_all(dapi_cent(:,1)<xmedline)'>manual_th];
[pos_cells_mat_sorted,XI] = sortrows(pos_cells_mat);
[uni_comb,ia] = unique(pos_cells_mat_sorted,'rows','first');
[~,ib] = unique(pos_cells_mat_sorted,'rows','last');
legend_str = repmat({''},length(uni_comb(:,1))-1,1);
for i=1:length(legend_str)
    if uni_comb(i+1,1)==1
        legend_str{i} = [legend_str{i},'&',chan1name];
    end
    if uni_comb(i+1,2)==1
        legend_str{i} = [legend_str{i},'&',chan2name];
    end
    if uni_comb(i+1,3)==1
        legend_str{i} = [legend_str{i},'&',chan4name];
    end
    legend_str{i} = legend_str{i}(2:end);
end
sum_combination = zeros(1,length(legend_str));
for i=2:length(legend_str)+1
    indtmp = XI(ia(i):ib(i));
    sum_combination(i-1) = ib(i)-ia(i)+1;
end
table1 = [[{'dapi'};{sum(dapi_cent(:,1)<xmedline & area_cell'>100)}],[legend_str';m2c(sum_combination)]];
saveCellFile(table1,[path_folder,filename(1:end-4),'_Left_counts_',date,'.txt']);

pos_cells_mat = [geneC1_dots_all(dapi_cent(:,1)>xmedline)'>manual_th,...
    geneC2_dots_all(dapi_cent(:,1)>xmedline)'>manual_th, geneC4_dots_all(dapi_cent(:,1)>xmedline)'>manual_th];
[pos_cells_mat_sorted,XI] = sortrows(pos_cells_mat);
[uni_comb,ia] = unique(pos_cells_mat_sorted,'rows','first');
[~,ib] = unique(pos_cells_mat_sorted,'rows','last');
legend_str = repmat({''},length(uni_comb(:,1))-1,1);
for i=1:length(legend_str)
    if uni_comb(i+1,1)==1
        legend_str{i} = [legend_str{i},'&',chan1name];
    end
    if uni_comb(i+1,2)==1
        legend_str{i} = [legend_str{i},'&',chan2name];
    end
    if uni_comb(i+1,3)==1
        legend_str{i} = [legend_str{i},'&',chan4name];
    end
    legend_str{i} = legend_str{i}(2:end);
end
sum_combination = zeros(1,length(legend_str));
for i=2:length(legend_str)+1
    indtmp = XI(ia(i):ib(i));
    sum_combination(i-1) = ib(i)-ia(i)+1;
end
table1 = [[{'dapi'};{sum(dapi_cent(:,1)>xmedline & area_cell'>100)}],[legend_str';m2c(sum_combination)]];
saveCellFile(table1,[path_folder,filename(1:end-4),'_Right_counts_',date,'.txt']);









