tic
clear all
close all


load /data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun20/afterLoading_analysis_dorsalhorn_all_fromdatabase_20-Jun-2016.mat
% % % % % % % % % % % % % 
geneuni = genes_dorsalhorn_all;

total_mol = sum(moldata_dorsalhorn_all);

uni_platevec = unique(chip_dorsalhorn_all);
samples_all = cellid_dorsalhorn_all;

% % % % % exclude low wells
goodwells = total_mol>1.5e3 & total_mol<70e3;% & ~strcmpi(annot_table_dorsalhorn_all.tissue,'hypothalamus') & ~strcmpi(annot_table_dorsalhorn_all.tissue,'striatum')
moldata_1 = moldata_dorsalhorn_all(:,goodwells);
chip_1 = chip_dorsalhorn_all(:,goodwells);
well_1 = well_dorsalhorn_all(:,goodwells);
age_1 = age_dorsalhorn_all(:,goodwells);
sex_1 = sex_dorsalhorn_all(:,goodwells);%% female=1, male =-1, mix=0;
diam_1 = diam_dorsalhorn_all(:,goodwells);
tissue_all_1 = tissue_dorsalhorn_all(goodwells);
strain_all_1 = strain_dorsalhorn_all(goodwells);
cellid_1 = cellid_dorsalhorn_all(goodwells);
total_mol_1 = total_mol(:,goodwells);
image_all_1 = image_all_dorsalhorn_all(goodwells);
moldata_spikes = moldata_dorsalhorn_all_ercc(:,goodwells);

tissue_uni = unique(tissue_all_1);

hbbgene = cellfun(@(x) ~isempty(strfind(x,'Hbb-')),geneuni(:,1)) | cellfun(@(x) ~isempty(strfind(x,'Hba-')),geneuni(:,1));
indgood = mean(moldata_1,2)>(50/length(chip_1)) & sum(moldata_1>0,2)<0.5*length(chip_1) & ~hbbgene;
moldata_1 = moldata_1(indgood,:);
geneuni_1 = geneuni(indgood,:);
% cellid_1 = cell(length(chip_1),1);
% for i=1:length(chip_1)
%     cellid_1{i} = [chip_1{i},'_',well_1{i}];
% end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%saveCellFile(table,['markertable_dorsalhorn_all',num2str(splitlev),'_',date,'.txt'])


% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

tic
splitlev = 5;
Nfeture = 100;
Nfeture1 = 500;
N_to_backspin = 10;
N_to_cells = 500;
mean_tresh = 0.01;
fdr_th = 0.3;
min_gr_cells = 5;
min_gr_genes = 20;
stop_th = [0.5,0.5];
flag_val_stip = 2;

% permute the cells order
% rng(1);
% cell_perm = randperm(length(cellid_1));

[dataout_sorted,genes_order,cells_order,genes_gr_level,cells_gr_level,genes_bor_level,cells_bor_level] =...
    backSpinSplit_v7(moldata_1,splitlev,Nfeture,Nfeture1,N_to_backspin,N_to_cells,mean_tresh,fdr_th,min_gr_cells,min_gr_genes,stop_th,flag_val_stip);



% [dataout_sorted,genes_order,cells_order,genes_gr_level,cells_gr_level,genes_bor_level,cells_bor_level] =...
%     backSpinSplit_v2(moldata_1(:,cell_perm),splitlev,Nfeture,Nfeture1,N_to_backspin,N_to_cells,mean_tresh,fdr_th,min_gr_cells,min_gr_genes);
toc


genes_sorted = geneuni_1(genes_order);
cells_sorted = cellid_1((cells_order));

save(['data_allcells_after_backspinv7_lev',num2str(splitlev),'_dorsalhorn_all_',date],'dataout_sorted', 'cells_sorted', 'genes_sorted'...
    ,'genes_bor_level','cells_bor_level','genes_gr_level','cells_gr_level','splitlev',....
    'Nfeture','N_to_backspin','N_to_cells','mean_tresh','fdr_th','min_gr_cells','min_gr_genes','stop_th','flag_val_stip');


% % load data_after_backspinv2_with_vm_lev5_dorsalhorn_all_23-Feb-2015.mat
sorted_data_cn = cent_norm(log2_center(dataout_sorted));

T_cells_tmp = cells_gr_level(:,splitlev+1)';
T_genes_tmp = genes_gr_level(:,splitlev+1)';




cells_bor_tmpgr = cells_bor_level;
genes_bor_tmpgr = genes_bor_level;

lev = splitlev;
cells_bor_tmpgr{lev} = [cells_bor_tmpgr{lev};length(T_cells_tmp)+1];
genes_bor_tmpgr{lev} = [genes_bor_tmpgr{lev};length(T_genes_tmp)+1];

figure;
set(gcf,'position',[100,100,450,750],'color','w')
axes('position',[0.1,0.02,0.88,0.85])
imagesc(sorted_data_cn,[prctile(sorted_data_cn(:),1),prctile(sorted_data_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
for jj=1:length(cells_bor_tmpgr{lev})
    if jj>1
        plot(cells_bor_tmpgr{lev}(jj)*[1,1]-0.5,[genes_bor_tmpgr{lev}(jj-1),genes_bor_tmpgr{lev}(jj)],'--','linewidth',linewid,'color',get_RGB(bor_color))
        plot([cells_bor_tmpgr{lev}(jj-1),cells_bor_tmpgr{lev}(jj)]-0.5,genes_bor_tmpgr{lev}(jj)*[1,1] ,'--','linewidth',linewid,'color',get_RGB(bor_color))
        plot(cells_bor_tmpgr{lev}(jj-1)*[1,1]-0.5,[genes_bor_tmpgr{lev}(jj-1),genes_bor_tmpgr{lev}(jj)]-0.5,'--','linewidth',linewid,'color',get_RGB(bor_color))
        plot([cells_bor_tmpgr{lev}(jj-1)-0.5,cells_bor_tmpgr{lev}(jj)],genes_bor_tmpgr{lev}(jj-1)*[1,1]-0.5 ,'--','linewidth',linewid,'color',get_RGB(bor_color))
    else
        plot(cells_bor_tmpgr{lev}(jj)*[1,1]-0.5,[1,genes_bor_tmpgr{lev}(jj)]-0.5,'--','linewidth',linewid,'color',get_RGB(bor_color))
        plot([1,cells_bor_tmpgr{lev}(jj)]-0.5,genes_bor_tmpgr{lev}(jj)*[1,1] -0.5,'--','linewidth',linewid,'color',get_RGB(bor_color))
        plot(1*[1,1]-0.5,[1,genes_bor_tmpgr{lev}(jj)]-0.5,'--','linewidth',linewid,'color',get_RGB(bor_color))
        plot([1,cells_bor_tmpgr{lev}(jj)]-0.5,1*[1,1] -0.5,'--','linewidth',linewid,'color',get_RGB(bor_color))
    end
end
set(gca,'xtick',[0:200:length(cells_sorted)],'ytick',[0:100:length(genes_sorted)], 'fontsize', 6)
colormap('summer');
freezeColors(gca);

[~,loc] = ismember(cells_sorted,cellid_1);
dataout_sorted_all = moldata_1(:,loc);


sorted_age_1 = age_1(loc);
sorted_sex_1 = sex_1(loc);
sorted_strain_1 = strain_all_1(loc);
sorted_diam_1 = diam_1(loc);
sorted_total_mol_1 = total_mol_1(loc);
sorted_chip_1 = chip_1(loc);
sorted_chip_1_uni = unique(sorted_chip_1);
sorted_chip_num = zeros(size(chip_1));
for i=1:length(sorted_chip_1_uni)
    sorted_chip_num(strcmpi(sorted_chip_1,sorted_chip_1_uni{i})) = i;
end

sorted_strain_1_uni = unique(sorted_strain_1);
sorted_strain_num = zeros(size(chip_1));
for i=1:length(sorted_strain_1_uni)
    sorted_strain_num(strcmpi(sorted_strain_1,sorted_strain_1_uni{i})) = i;
end

sorted_tissue_1 = tissue_all_1(loc);
sorted_tissue_num = zeros(size(tissue_all_1));
sorted_tissue_nummat = zeros(length(tissue_uni),length(tissue_all_1));
for i=1:length(tissue_uni)
    sorted_tissue_num(strcmpi(sorted_tissue_1,tissue_uni{i})) = i;
    sorted_tissue_nummat(i, strcmpi(sorted_tissue_1,tissue_uni{i})) = 1;
end

sorted_age_1(isnan(sorted_age_1)) = mean(sorted_age_1(~isnan(sorted_age_1)));

for i=1:length(tissue_uni)
    axes('position',[0.1,0.87+(i-1)*0.005,0.88,0.005]);
    imagesc(-sorted_tissue_nummat(i,:));
    set(gca,'xtick',[],'ytick',[]);
    text(-0.5,1,tissue_uni{i},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',6);
    colormap('gray');
    freezeColors(gca);
end

% axes('position',[0.1,0.88,0.85,0.01])
% imagesc(cent_norm(T_cells_tmp))
% set(gca,'xtick',[],'ytick',1,'yticklabel','clusters')
axes('position',[0.1,0.96,0.88,0.007])
imagesc(cent_norm(sorted_age_1))
set(gca,'xtick',[],'ytick',1,'yticklabel','age')
colormap('jet');freezeColors(gca);
% axes('position',[0.1,0.94,0.88,0.01])
% imagesc(cent_norm(sorted_chip_num))
% set(gca,'xtick',[],'ytick',1,'yticklabel','chip')
axes('position',[0.1,0.967,0.88,0.007])
imagesc(cent_norm(sorted_sex_1))
set(gca,'xtick',[],'ytick',1,'yticklabel','sex')
colormap('jet');freezeColors(gca);
axes('position',[0.1,0.974,0.88,0.007])
imagesc(cent_norm(sorted_diam_1))
set(gca,'xtick',[],'ytick',1,'yticklabel','diam');
colormap('jet');freezeColors(gca);
% axes('position',[0.1,0.96,0.88,0.01])
% imagesc(cent_norm(sorted_tissue_num))
% set(gca,'xtick',[],'ytick',1,'yticklabel','tissue')

axes('position',[0.1,0.981,0.88,0.01])
imagesc(cent_norm(log2(sorted_total_mol_1)))
set(gca,'xtick',[],'ytick',1,'yticklabel','tot. mol.')
colormap('jet');freezeColors(gca);

axes('position',[0.1,0.92,0.88,0.01])
imagesc(cent_norm(log2(sorted_strain_num)))
set(gca,'xtick',[],'ytick',1,'yticklabel','strain')
colormap('colorcube');freezeColors(gca);


eval(['export_fig -r2000 allcells_lev',num2str(splitlev),'_',date,'.pdf ']);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % %
T_cells_tmp_uni = unique(T_cells_tmp);
gr_center = zeros(length(T_cells_tmp_uni),1);
% % % % % % % % % % % % % % % % 

meanpergene = mean(dataout_sorted_all,2);
molenrich_mat = zeros(length(dataout_sorted_all(:,1)),length(T_cells_tmp_uni));
meangr_mat = zeros(length(dataout_sorted_all(:,1)),length(T_cells_tmp_uni));
meangrpos_mat = zeros(length(dataout_sorted_all(:,1)),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
    molenrich_mat(:,jjj) = mean(dataout_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./meanpergene;
    meangrpos_mat(:,jjj) = mean(dataout_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj))>0,2);
end
% molenrich_mat(meangrpos_mat<0.1) = 0;
meangrpos_mat(meangrpos_mat<0.2) = 0;
molenrich_mat(meanpergene==0 | isnan(meanpergene),:) = 0;
molenrich_mat = molenrich_mat.*meangrpos_mat.^1;
[~,grord] = sort(gr_center);
molenrich_mat = molenrich_mat(:,grord);
[~,xi] = sort(molenrich_mat,'descend');


ind_gr_tmp_mark = xi(1:2,:);
ind_gr_tmp_mark = flipud(ind_gr_tmp_mark(:));
rmv = false(size(ind_gr_tmp_mark));
for i=2:length(ind_gr_tmp_mark)
    if ismember(ind_gr_tmp_mark(i),ind_gr_tmp_mark(1:i-1))
        rmv(i) = true;
    end
end
ind_gr_tmp_mark(rmv) = [];
gr_tmp_mark = geneuni_1(ind_gr_tmp_mark);
gr_tmp_mark = (gr_tmp_mark(:));

figure;
set(gcf,'position',[100,100,450,750],'color','w')

barhight = 0.89/length(ind_gr_tmp_mark);
for jj=1:length(ind_gr_tmp_mark)
    disp([num2str(jj),','])
    maxval = prctile(dataout_sorted_all(ind_gr_tmp_mark(jj),:),99.99);%0.9*max(data_tmpgr(ind_gr_tmp_mark(jj),:));
    n = length(T_cells_tmp);
    axes('position',[0.12,0.01+barhight*(jj-1),0.88,barhight])
    %     if mod(jj,2)==1;
    %         bar(n/2,maxval,'facecolor',get_RGB('grey81'),'edgecolor','none','BarWidth',n);hold on;axis tight;axis off;
    %     end
    if mod(jj,2)==1;
        h = bar(dataout_sorted_all(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',1);hold on;axis tight;axis off;
    else
        h = bar(dataout_sorted_all(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',1,'facecolor',[200,0,0]/256);hold on;axis tight;axis off;
    end
    for jjj=1:length(cells_bor_tmpgr{lev})
        plot(cells_bor_tmpgr{lev}(jjj)*[1,1]-1,[0,maxval],'--','linewidth',0.5,'color',[200,200,200]/256); hold on;
    end
    
    set(gca,'xlim',[0,n],'ylim',[0,maxval], 'fontsize', 6)
    text(-0.5,maxval/2,gr_tmp_mark{jj},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',7)
end


eval(['export_fig barplots_allcells_top2_sum_backspinv2_lev',num2str(splitlev),'_',date,'.pdf'])

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
sumpergene = sum(dataout_sorted_all,2);
meanpergene = mean(dataout_sorted_all,2);%
molenrich_mat = zeros(length(dataout_sorted_all(:,1)),length(T_cells_tmp_uni));
meangr_mat = zeros(length(dataout_sorted_all(:,1)),length(T_cells_tmp_uni));
meangrpos_mat = zeros(length(dataout_sorted_all(:,1)),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
    molenrich_mat(:,jjj) = mean(dataout_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./meanpergene;
    meangrpos_mat(:,jjj) = mean(dataout_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj))>0,2);
end
[~,grord] = sort(gr_center);
meangrpos_mat = meangrpos_mat(:,grord);
molenrich_mat = molenrich_mat(:,grord);
[~,xi0] = sort(molenrich_mat.*meangrpos_mat.^0,'descend');
[~,xi0p5] = sort(molenrich_mat.*meangrpos_mat.^0.5,'descend');
[~,xi1] = sort(molenrich_mat.*meangrpos_mat.^1,'descend');


ind_gr_tmp_mark = [xi0(1:5,:);xi0p5(1:5,:);xi1(1:5,:)];
ind_gr_tmp_mark = flipud(ind_gr_tmp_mark(:));
rmv = false(size(ind_gr_tmp_mark));
for i=2:length(ind_gr_tmp_mark)
    if ismember(ind_gr_tmp_mark(i),ind_gr_tmp_mark(1:i-1))
        rmv(i) = true;
    end
end
ind_gr_tmp_mark(rmv) = [];
gr_tmp_mark = geneuni_1(ind_gr_tmp_mark);
gr_tmp_mark = (gr_tmp_mark(:));

datamarkers = dataout_sorted_all(ind_gr_tmp_mark,:);
datamarkers_cn = cent_norm(log2(datamarkers+1));


figure;
set(gcf,'position',[100,100,450,750],'color','w')
axes('position',[0.1,0.02,0.88,0.85])
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
for jj=1:length(cells_bor_tmpgr{lev})
    plot(cells_bor_tmpgr{lev}(jj)*[1,1]-0.5,[1,length(gr_tmp_mark)],'-','linewidth',linewid,'color',get_RGB(bor_color))
end
set(gca,'xtick',[0:200:length(cells_sorted)],'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 2)
colormap('summer');
freezeColors(gca);


for i=1:length(tissue_uni)
    axes('position',[0.1,0.87+(i-1)*0.006,0.88,0.006]);
    imagesc(-sorted_tissue_nummat(i,:));
    set(gca,'xtick',[],'ytick',[], 'fontsize', 6);
    text(-0.5,1,tissue_uni{i},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',6);
    colormap('gray');
    freezeColors(gca);
end

eval(['export_fig -r2000 heatmap_for_allcells_markertable_backspinv7_lev',num2str(splitlev),'_',date,'.pdf'])

table = m2c([T_cells_tmp;datamarkers]');
table = [[{'','strain','tissue','cluster'},gr_tmp_mark'];[cells_sorted',sorted_strain_1',sorted_tissue_1',table]];

saveCellFile(table,['markertable_allcells_backspinv7_lev',num2str(splitlev),'_',date,'.txt'])





