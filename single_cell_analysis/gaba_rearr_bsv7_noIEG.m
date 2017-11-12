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
goodwells = total_mol>3e3 & total_mol<70e3;% & ~strcmpi(annot_table_dorsalhorn_all.tissue,'hypothalamus') & ~strcmpi(annot_table_dorsalhorn_all.tissue,'striatum')
moldata_1 = moldata_dorsalhorn_all(:,goodwells);
chip_1 = chip_dorsalhorn_all(:,goodwells);
well_1 = well_dorsalhorn_all(:,goodwells);
age_1 = age_dorsalhorn_all(:,goodwells);
sex_1 = sex_dorsalhorn_all(:,goodwells);%% female=1, male =-1, mix=0;
diam_1 = diam_dorsalhorn_all(:,goodwells);
tissue_all_1 = tissue_dorsalhorn_all(goodwells);
total_mol_1 = total_mol(:,goodwells);
image_all_1 = image_all_dorsalhorn_all(goodwells);
green_1 = green_dorsalhorn_all(:,goodwells);
strain_1 = strain_dorsalhorn_all(:,goodwells);
cellid_1 = cellid_dorsalhorn_all(:,goodwells);
donorid_1 = donorid_dorsalhorn_all(:,goodwells);


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




% % % % % % % % % % % % % % % % % % % % % % %
%Extracting neurons


markertable = loadCellFile_turbo('lev1_markertable_allcells_backspinv7_lev5_20-Jun-2016.txt',1);

markertable = markertable(2:end,1:5);

neuronsind = find(strcmpi(markertable(:,5),'neuron'));

cellid_neurons = markertable(neuronsind,1);


[~,loc_neurons] = ismember(cellid_neurons, cellid_1);
loc_neurons(loc_neurons==0) = [];

%loc_nonneurons = setdiff(cellid_1, cellid_neurons
loc_nonneurons = setdiff([1:length(cellid_1)], loc_neurons);




% % % % % % % % % % % % % % %
% t-test to exclude non-neuron genes



[~,p] = ttest2(log2(moldata_1(:,loc_neurons)+1)',log2(moldata_1(:,loc_nonneurons)+1)','tail','right');

%figure;plot(sort(p))

neurons_genes_ind = intersect(setdiff([1:length(geneuni_1)],fdr_proc(1-p,0.1)),find(sum(moldata_1(:,loc_neurons),2)>10));


%find(strcmpi(geneuni_1(neurons_genes_ind),'aqp4'))
%find(strcmpi(geneuni_1(neurons_genes_ind),'sst'))





% % % % % % % % % % % % % % %
% exclude IEG

geneuni_1 = geneuni_1(neurons_genes_ind);


ie_genes = loadCellFile('IEG_mod.txt');

[~,loc] = ismember(ie_genes,geneuni_1);
ieg_ind = loc(loc>0);


neurons_genes_ind(ieg_ind) = [];

geneuni_1(ieg_ind) = [];




% % % % % % % % % % % % % % % % % % % % % % %
%Extracting gaba


neuronsmarkertable = loadCellFile_turbo('classif_hannah_markertable_neurons_backspinv2_lev7_21-Jun-2016.txt',1);

neuronsmarkertable = neuronsmarkertable(2:end,1:4);

gabaind = find(strcmpi(neuronsmarkertable(:,4),'gaba'));

cellid_gaba = neuronsmarkertable(gabaind,1);


[~,loc_gaba] = ismember(cellid_gaba, cellid_1);
loc_gaba(loc_gaba==0) = [];



% % % % % % % % % % % % % % % %
% 

moldata_1 = moldata_1(neurons_genes_ind,loc_gaba);
chip_1 = chip_1(loc_gaba);
well_1 = well_1(loc_gaba);
age_1 = age_1(loc_gaba);
sex_1 = sex_1(loc_gaba);
diam_1 = diam_1(loc_gaba);
tissue_all_1 = tissue_all_1(loc_gaba);
total_mol_1 = total_mol_1(loc_gaba);
image_all_1 = image_all_1(loc_gaba);
moldata_spikes = moldata_spikes(:,loc_gaba);
green_1 = green_1(loc_gaba);
strain_1 = strain_1(loc_gaba);
cellid_1 = cellid_1(loc_gaba);
donorid_1 = donorid_1(loc_gaba);




% redefine boarders after classification

markertable = loadCellFile_turbo('/data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/classif_gaba_noIEG_bsv7_lev5_22-Jul-2016.txt',1);
markertable = markertable(2:end,1:6);

lev2_num = cellfun(@str2double, markertable(:,4));
[~,xi] = sort(lev2_num);

markertable = markertable(xi,:);
lev2_num = lev2_num(xi,:);

[~,loc] = ismember(markertable(:,1), cellid_1);
dataout_sorted_all = moldata_1(:,loc);
sorted_strain_1 = strain_1(loc);

sorted_strain_1_uni = unique(sorted_strain_1);
sorted_strain_num = zeros(size(chip_1));
for i=1:length(sorted_strain_1_uni)
    sorted_strain_num(strcmpi(sorted_strain_1,sorted_strain_1_uni{i})) = i;
end



cells_bor_tmpgr = find(diff(lev2_num)==1)+1;

T_cells_tmp = lev2_num;

cells_bor_tmpgr = [cells_bor_tmpgr;length(T_cells_tmp)+1];


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % %


T_cells_tmp_uni = unique(T_cells_tmp);
gr_center = zeros(length(T_cells_tmp_uni),1);

% % % % % % % % % % % % % % %
% Marker genes barplots



sumpergene = sum(dataout_sorted_all,2);
meanpergene = mean(dataout_sorted_all,2);%mean(log2(dataout_sorted_all+1),2);
molenrich_mat = zeros(length(dataout_sorted_all(:,1)),length(T_cells_tmp_uni));
meangr_mat = zeros(length(dataout_sorted_all(:,1)),length(T_cells_tmp_uni));
meangrpos_mat = zeros(length(dataout_sorted_all(:,1)),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
%     molenrich_mat(:,jjj) = sum(dataout_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./sumpergene*length(T_cells_tmp)/sum(T_cells_tmp==T_cells_tmp_uni(jjj));
    molenrich_mat(:,jjj) = mean(dataout_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./meanpergene;
%     meangr_mat(:,jjj) = mean(log2(dataout_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj))+1),2);
    meangrpos_mat(:,jjj) = mean(dataout_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj))>0,2);
end
% molenrich_mat(meangrpos_mat<0.1) = 0;
meangrpos_mat(meangrpos_mat<0.2) = 0;
molenrich_mat(meanpergene==0,:) = 0;
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
set(gcf,'position',[100,100,550,1000],'color','w')
barhight = 0.95/length(ind_gr_tmp_mark);
for jj=1:length(ind_gr_tmp_mark)
    disp([num2str(jj),','])
    maxval = prctile(dataout_sorted_all(ind_gr_tmp_mark(jj),:),99.99);%0.9*max(data_tmpgr(ind_gr_tmp_mark(jj),:));
    n = length(T_cells_tmp);
    axes('position',[0.12,0.01+barhight*(jj-1),0.88,barhight])
    %     if mod(jj,2)==1;
    %         bar(n/2,maxval,'facecolor',get_RGB('grey81'),'edgecolor','none','BarWidth',n);hold on;axis tight;axis off;
    %     end
    if mod(jj,2)==1;
        h = bar(dataout_sorted_all(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',3);hold on;axis tight;axis off;
    else
        h = bar(dataout_sorted_all(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',3,'facecolor',[200,0,0]/256);hold on;axis tight;axis off;
    end
    for jjj=1:length(cells_bor_tmpgr)
        plot(cells_bor_tmpgr(jjj)*[1,1]-1,[0,maxval],'--','linewidth',0.5,'color',[200,200,200]/256); hold on;
    end
    
    set(gca,'xlim',[0,n],'ylim',[0,maxval], 'fontsize', 6)
    text(-0.5,maxval/2,gr_tmp_mark{jj},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',7)
    if jj==1
        axis on; box off;
     set(gca,'xtick',round(gr_center(grord)),'xticklabel',[1:length(gr_center)],'ytick',[]);
    end
end

axes('position',[0.12,0.98,0.88,0.01])
imagesc(cent_norm(log2(sorted_strain_num)))
set(gca,'xtick',[],'ytick',1,'yticklabel','strain')
colormap('colorcube');freezeColors(gca);


% eval(['export_fig barplots_gaba_rearr_noIEG_top2_sum_bsv7_',date,'.pdf'])


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% barplots with selected genes


gr_tmp_mark = loadCellFile('new_list_gaba.txt');
[~,ind_gr_tmp_mark] = ismember(gr_tmp_mark,geneuni_1);
gr_tmp_mark(ind_gr_tmp_mark==0) = [];
[~,ind_gr_tmp_mark] = ismember(gr_tmp_mark,geneuni_1);

figure;
set(gcf,'position',[100,100,450,750],'color','w')

barhight = 0.89/length(ind_gr_tmp_mark);
for jj=1:length(ind_gr_tmp_mark)
   disp([num2str(jj),','])
   maxval = prctile(dataout_sorted_all(ind_gr_tmp_mark(jj),:),99.99);%0.9*max(data_tmpgr(ind_gr_tmp_mark(jj),:));
   n = length(T_cells_tmp);
   axes('position',[0.12,0.01+barhight*(jj-1),0.8,barhight])
   %     if mod(jj,2)==1;
   % bar(n/2,maxval,'facecolor',get_RGB('grey81'),'edgecolor','none','BarWidth',n);hold on;axis tight;axis off;
   %     end
   if mod(jj,2)==1;
       h = bar(dataout_sorted_all(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',3);hold on;axis tight;axis off;
   else
       h = bar(dataout_sorted_all(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',3,'facecolor',[200,0,0]/256);hold on;axis tight;axis off;
   end
   for jjj=1:length(cells_bor_tmpgr)
plot(cells_bor_tmpgr(jjj)*[1,1]-1,[0,maxval],'--','linewidth',0.5,'color',[200,200,200]/256); hold on;
   end

   set(gca,'xlim',[0,n],'ylim',[0,maxval], 'fontsize', 6)
text(-0.5,maxval/2,gr_tmp_mark{jj},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',7)
end


% eval(['export_fig barplots_gabanoIEGmarkers_rearr_',date,'.pdf'])
% 





% % % % % % % % % % % % % %
% Neurotransmitters & Receptors



gr_tmp_mark = loadCellFile('/data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/neurotransm_receptors.txt');
[~,ind_gr_tmp_mark] = ismember(gr_tmp_mark,geneuni_1);
gr_tmp_mark(ind_gr_tmp_mark==0) = [];
[~,ind_gr_tmp_mark] = ismember(gr_tmp_mark,geneuni_1);

figure;
set(gcf,'position',[100,100,450,750],'color','w')

barhight = 0.89/length(ind_gr_tmp_mark);
for jj=1:length(ind_gr_tmp_mark)
   disp([num2str(jj),','])
   maxval = prctile(dataout_sorted_all(ind_gr_tmp_mark(jj),:),99.99);%0.9*max(data_tmpgr(ind_gr_tmp_mark(jj),:));
   n = length(T_cells_tmp);
   axes('position',[0.12,0.01+barhight*(jj-1),0.8,barhight])
   %     if mod(jj,2)==1;
   % bar(n/2,maxval,'facecolor',get_RGB('grey81'),'edgecolor','none','BarWidth',n);hold on;axis tight;axis off;
   %     end
   if mod(jj,2)==1;
       h = bar(dataout_sorted_all(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',3);hold on;axis tight;axis off;
   else
       h = bar(dataout_sorted_all(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',3,'facecolor',[200,0,0]/256);hold on;axis tight;axis off;
   end
   for jjj=1:length(cells_bor_tmpgr)
plot(cells_bor_tmpgr(jjj)*[1,1]-1,[0,maxval],'--','linewidth',0.5,'color',[200,200,200]/256); hold on;
   end

   set(gca,'xlim',[0,n],'ylim',[0,maxval], 'fontsize', 6)
text(-0.5,maxval/2,gr_tmp_mark{jj},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',7)
end


eval(['export_fig barplots_gaba_rearr_neurotransm_noIEG_bsv7_',date,'.pdf'])



% % % % % % % % % % % % % %
% Peptides & Receptors


gr_tmp_mark = loadCellFile('/data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/peptides_receptors.txt');
[~,ind_gr_tmp_mark] = ismember(gr_tmp_mark,geneuni_1);
gr_tmp_mark(ind_gr_tmp_mark==0) = [];
[~,ind_gr_tmp_mark] = ismember(gr_tmp_mark,geneuni_1);

figure;
set(gcf,'position',[100,100,450,750],'color','w')

barhight = 0.89/length(ind_gr_tmp_mark);
for jj=1:length(ind_gr_tmp_mark)
   disp([num2str(jj),','])
   maxval = prctile(dataout_sorted_all(ind_gr_tmp_mark(jj),:),99.99);%0.9*max(data_tmpgr(ind_gr_tmp_mark(jj),:));
   n = length(T_cells_tmp);
   axes('position',[0.12,0.01+barhight*(jj-1),0.8,barhight])
   %     if mod(jj,2)==1;
   % bar(n/2,maxval,'facecolor',get_RGB('grey81'),'edgecolor','none','BarWidth',n);hold on;axis tight;axis off;
   %     end
   if mod(jj,2)==1;
       h = bar(dataout_sorted_all(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',3);hold on;axis tight;axis off;
   else
       h = bar(dataout_sorted_all(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',3,'facecolor',[200,0,0]/256);hold on;axis tight;axis off;
   end
   for jjj=1:length(cells_bor_tmpgr)
plot(cells_bor_tmpgr(jjj)*[1,1]-1,[0,maxval],'--','linewidth',0.5,'color',[200,200,200]/256); hold on;
   end

   set(gca,'xlim',[0,n],'ylim',[0,maxval], 'fontsize', 6)
text(-0.5,maxval/2,gr_tmp_mark{jj},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',7)
end



eval(['export_fig barplots_gaba_rearr_peptides_noIEG_bsv7_lev',date,'.pdf'])





% % % % % % % % % % % % % % % %


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

meangrpos_mat(meangrpos_mat<0.2) = 0;
molenrich_mat(meanpergene==0,:) = 0;

% meangrpos_mat(meangrpos_mat<0.2) = 0;
meangrpos_mat = meangrpos_mat(:,grord);
molenrich_mat = molenrich_mat(:,grord);
[~,xi0] = sort(molenrich_mat.*meangrpos_mat.^0.001,'descend');
[~,xi0p5] = sort(molenrich_mat.*meangrpos_mat.^0.5,'descend');
[~,xi1] = sort(molenrich_mat.*meangrpos_mat.^1,'descend');


ind_gr_tmp_mark = [xi0(1:10,:);xi0p5(1:10,:);xi1(1:10,:)];
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
set(gcf,'position',[100,100,550,1000],'color','w')
axes('position',[0.1,0.02,0.88,0.95])
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
for jj=1:length(cells_bor_tmpgr)
    plot(cells_bor_tmpgr(jj)*[1,1]-0.5,[1,length(gr_tmp_mark)],'-','linewidth',linewid,'color',get_RGB(bor_color))
end
set(gca,'xtick',[0:200:length(cellid_1)],'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 2)
colormap('summer');
freezeColors(gca);
set(gca,'xtick',round(gr_center(grord)),'xticklabel',[1:length(gr_center)]);

axes('position',[0.1,0.98,0.88,0.01])
imagesc(cent_norm(log2(sorted_strain_num)))
set(gca,'xtick',[],'ytick',1,'yticklabel','strain')
colormap('colorcube');freezeColors(gca);



% eval(['export_fig heatmap_gaba_rearr_noIEG_markertable_bsv7_',date,'.pdf'])




% % %save2pdf(['heatmap_for_gaba_markertable_backspinv2_lev',num2str(splitlev),'_',date,'.pdf'],gcf,300) 
% 
% 
% 
% % % % % % % % 
% % save marker table
% 
% table = m2c([T_cells_tmp;datamarkers]');
% table = [[{'','tissue','cluster'},gr_tmp_mark'];[cells_sorted',sorted_tissue_1',table]];
% 
% 
% saveCellFile(table,['markertable_gaba_noIEG_backspinv2_lev',num2str(splitlev),'_',date,'.txt'])
% 
% 
% 
% % 
% % 
% % 
% 
