
tic
clear all
close all


load /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun20/afterLoading_analysis_dorsalhorn_all_fromdatabase_20-Jun-2016.mat
% % % % % % % % % % % % % 
geneuni = genes_dorsalhorn_all;

total_mol = sum(moldata_dorsalhorn_all);

uni_platevec = unique(chip_dorsalhorn_all);
samples_all = cellid_dorsalhorn_all;

% % % % % exclude low wells
goodwells = total_mol>3.0e3 & total_mol<70e3;% & ~strcmpi(annot_table_dorsalhorn_all.tissue,'hypothalamus') & ~strcmpi(annot_table_dorsalhorn_all.tissue,'striatum')
moldata_1 = moldata_dorsalhorn_all(:,goodwells);
chip_1 = chip_dorsalhorn_all(:,goodwells);
well_1 = well_dorsalhorn_all(:,goodwells);
age_1 = age_dorsalhorn_all(:,goodwells);
sex_1 = sex_dorsalhorn_all(:,goodwells);%% female=1, male =-1, mix=0;
diam_1 = diam_dorsalhorn_all(:,goodwells);
tissue_all_1 = tissue_dorsalhorn_all(goodwells);
total_mol_1 = total_mol(:,goodwells);
image_all_1 = image_all_dorsalhorn_all(goodwells);
moldata_spikes = moldata_dorsalhorn_all_ercc(:,goodwells);

tissue_uni = unique(tissue_all_1);

% hbbgene = cellfun(@(x) ~isempty(strfind(x,'Hbb-')),geneuni(:,1)) | cellfun(@(x) ~isempty(strfind(x,'Hba-')),geneuni(:,1));
indgood = mean(moldata_1,2)>(10/length(chip_1)) & sum(moldata_1>0,2)<1*length(chip_1);
moldata_1 = moldata_1(indgood,:);
geneuni_1 = geneuni(indgood,:);
cellid_1 = cell(length(chip_1),1);
for i=1:length(chip_1)
    cellid_1{i} = [chip_1{i},'_',well_1{i}];
end


% % % % % % % % % % % % % % % % % % % % % % %
%Extracting neurons
markertable = loadCellFile_turbo('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/lev1_markertable_allcells_backspinv7_lev5_20-Jun-2016.txt',1);
markertable = markertable(2:end,1:5);
neuronsind = find(strcmpi(markertable(:,5),'neuron'));
cellid_neurons = markertable(neuronsind,1);
[~,loc_neurons] = ismember(cellid_neurons, cellid_1);
loc_neurons(loc_neurons==0) = [];

%loc_nonneurons = setdiff(cellid_1, cellid_neurons
loc_nonneurons = setdiff([1:length(cellid_1)], loc_neurons);



% % % % % % % % % % % % % % % %
% % t-test to exclude non-neuron genes
% [~,p] = ttest2(log2(moldata_1(:,loc_neurons)+1)',log2(moldata_1(:,loc_nonneurons)+1)','tail','right');
% %figure;plot(sort(p))
% neurons_genes_ind = intersect(setdiff([1:length(geneuni_1)],fdr_proc(1-p,0.5)),find(sum(moldata_1(:,loc_neurons),2)>10));
% %find(strcmpi(geneuni_1(neurons_genes_ind),'aqp4'))
% %find(strcmpi(geneuni_1(neurons_genes_ind),'sst'))
% % % % % % % % % % % % % % % %
% geneuni_1 = geneuni_1(neurons_genes_ind);


moldata_1 = moldata_1(:,loc_neurons);
chip_1 = chip_1(loc_neurons);
well_1 = well_1(loc_neurons);
age_1 = age_1(loc_neurons);
sex_1 = sex_1(loc_neurons);
diam_1 = diam_1(loc_neurons);
tissue_all_1 = tissue_all_1(loc_neurons);
total_mol_1 = total_mol_1(loc_neurons);
image_all_1 = image_all_1(loc_neurons);
moldata_spikes = moldata_spikes(:,loc_neurons);
cellid_1 = cell(length(chip_1),1);
for i=1:length(chip_1)
    cellid_1{i} = [chip_1{i},'_',well_1{i}];
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% load /data2/C1_stuff/DentGyr/data_after_backspinv2_nonneurons_lev4_DG_23-Dec-2015.mat
lev2_annot = loadCellFile_turbo('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/classif_all_noIEG_bsv7_lev5_22-Jul-2016_excl.txt',1);
lev2_annot = lev2_annot(2:end,1:7);
lev2_annot(strcmpi(lev2_annot(:,7),'Glut_exclude') | strcmpi(lev2_annot(:,7),'GABA_exclude'),:) = [];

grlev2 = lev2_annot(:,7);
grlev2_uni = unique(grlev2);

[~,loc] = ismember(lev2_annot(:,1) ,cellid_1);

moldata_1 = moldata_1(:,loc);
chip_1 = chip_1(loc);
well_1 = well_1(loc);
age_1 = age_1(loc);
sex_1 = sex_1(loc);
diam_1 = diam_1(loc);
tissue_all_1 = tissue_all_1(loc);
total_mol_1 = total_mol_1(loc);
image_all_1 = image_all_1(loc);
moldata_spikes = moldata_spikes(:,loc);
cellid_1 = cellid_1(loc);

dend_order = loadCellFile('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/lev2_dendrogram_order_nov25_2016.txt');
dend_order = regexprep(dend_order, ' ','');
dend_order = flipud(dend_order);
grlev2_uni = dend_order;

indgood = [1:length(moldata_1(1,:))];
indgood2 = indgood;
% indgood2 = find(sorted_total_mol_1>2500);
T_cells_tmp = zeros(length(grlev2),1);
for i=1:length(grlev2_uni)
    T_cells_tmp( strcmpi(grlev2, grlev2_uni{i}) ) = i;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % %
[T_cells_tmp,xi] = sort(T_cells_tmp);
xi(T_cells_tmp==0) = [];
T_cells_tmp(T_cells_tmp==0) = [];
T_cells_tmp_uni = unique(T_cells_tmp);

moldata_1 = moldata_1(:,xi);
chip_1 = chip_1(xi);
well_1 = well_1(xi);
age_1 = age_1(xi);
sex_1 = sex_1(xi);
diam_1 = diam_1(xi);
tissue_all_1 = tissue_all_1(xi);
total_mol_1 = total_mol_1(xi);
image_all_1 = image_all_1(xi);
moldata_spikes = moldata_spikes(:,xi);
cellid_1 = cellid_1(xi);

grlev2 = grlev2(xi);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
fname = dir('*png');
fname = {fname.name}';
fname = regexprep(fname,'_level2_group.png','');
genesleft = setdiff(geneuni_1(:,1),fname);
[~,loc] = ismember(genesleft,geneuni_1);
moldata_1 = moldata_1(loc,:);
geneuni_1 = geneuni_1(loc,1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
moldata_validcells = moldata_1;
dend_order_names = regexprep(dend_order,'_','-');
% Tcelluni = unique(lev2_gr_num_treeorder);
space_gr = 20;
% colorvec = distinguishable_colors(length(T_cells_tmp_uni));
% colorvec = distinguishable_colors(15);
% colorvec = [colorvec;colorvec];
colors_glut = [222 220 0
252 234 16
149 193 31
0 141 54
29 113 184
54 169 225
75 81 140
190 22 34
233 78 27
249 178 51
243 146 0
163 25 91
242 87 202
246 158 220
175 104 156]/255; 
color_gaba = [252 234 16
222 220 0
149 193 31
0 141 54
54 169 225
29 113 184
75 81 140
249 178 51
243 146 0
163 25 91
230 51 42
190 22 34
233 78 27
202 158 103
147 96 55]/255;
colorvec = [colors_glut;color_gaba];
for i=1:length(T_cells_tmp_uni)
    ind = find(T_cells_tmp==i);
    gr_cent(i) = mean(ind + (i-1)*space_gr);
end

tree_order_label = {'Glut1','Glut2','Glut3','Glut4','Glut5','Glut6','Glut7','Glut8','Glut9','Glut10','Glut11','Glut12','Glut13','Glut14','Glut15',.....
    'Gaba1','Gaba2','Gaba3','Gaba4','Gaba5','Gaba6','Gaba7','Gaba8','Gaba9','Gaba10','Gaba11','Gaba12','Gaba13','Gaba14','Gaba15'};

yax_lim = [0,space_gr*length(T_cells_tmp_uni) + length(moldata_validcells(1,:))];


hf = figure('color','w','position',[100,100,1200,300],'visible','on');
hold on;
totmol = sum(moldata_validcells);
maxval = 0.8*max(totmol);
totmol_mean = zeros(size(T_cells_tmp_uni));
totmol_std = zeros(size(T_cells_tmp_uni));
cellnum = zeros(size(T_cells_tmp_uni));

for i=1:length(T_cells_tmp_uni)
    ind = find(T_cells_tmp==i);
    cellnum(i) = length(ind);
    totmol_mean(i) = mean(totmol(ind));
    totmol_std(i) = std(totmol(ind));
    plot((ind(1) + (i-1)*space_gr-1)*[1,1],[0,maxval], '--','color',220/255*[1,1,1]);
    plot((ind(end) + (i-1)*space_gr+0.5)*[1,1],[0,maxval], '--','color',220/255*[1,1,1]);
    h = bar(ind + (i-1)*space_gr,totmol(ind));
    set(h,'facecolor',colorvec(i,:),'edgecolor',colorvec(i,:),'barwidth',1);
end
heb = errorbar(gr_cent, totmol_mean, totmol_std, '.');
set(heb,'linewidth',3,'capsize', 10);   
plot(gr_cent, totmol_mean, 'sk');
set(gca,'position',[0.05,0.3,0.9,0.6])
set(gca,'xtick',gr_cent,'xticklabel',tree_order_label,'fontsize',14,'ylim',[0,maxval],'xlim',yax_lim,'xdir','normal','XTickLabelRotation',45)
ylabel('mol/cell','fontsize',14);
title('Total mol.','fontsize',14);
set(gca,'xlim',yax_lim);

eval(['export_fig barplot_totmol_percell_',date,'.pdf'])


hf = figure('color','w','position',[100,100,1200,300],'visible','on');
hold on;
totgene = sum(moldata_validcells>0);
maxval = 0.8*max(totgene);
totgene_mean = zeros(size(T_cells_tmp_uni));
totgene_std = zeros(size(T_cells_tmp_uni));
cellnum = zeros(size(T_cells_tmp_uni));

for i=1:length(T_cells_tmp_uni)
    ind = find(T_cells_tmp==i);
    cellnum(i) = length(ind);
    totgene_mean(i) = mean(totgene(ind));
    totgene_std(i) = std(totgene(ind));
    plot((ind(1) + (i-1)*space_gr-1)*[1,1],[0,maxval], '--','color',220/255*[1,1,1]);
    plot((ind(end) + (i-1)*space_gr+0.5)*[1,1],[0,maxval], '--','color',220/255*[1,1,1]);
    h = bar(ind + (i-1)*space_gr,totgene(ind));
    set(h,'facecolor',colorvec(i,:),'edgecolor',colorvec(i,:),'barwidth',1);
end
heb = errorbar(gr_cent, totgene_mean, totgene_std, '.');
set(heb,'linewidth',3,'capsize', 10);   
plot(gr_cent, totgene_mean, 'sk');
set(gca,'position',[0.05,0.3,0.9,0.6])
set(gca,'xtick',gr_cent,'xticklabel',tree_order_label,'fontsize',14,'ylim',[0,maxval],'xlim',yax_lim,'xdir','normal','XTickLabelRotation',45)
ylabel('genes/cell','fontsize',14);
title('Total genes.','fontsize',14);
set(gca,'xlim',yax_lim);

eval(['export_fig barplot_totgenes_percell_',date,'.pdf'])


table1 = [tree_order_label;m2c([cellnum';totmol_mean';totmol_std';totgene_mean';totgene_std'])];
table1 = [{'clustername';'#cells';'totmol-mean';'totmol-std';'totgenes-mean';'totgenes-std'},table1];

saveCellFile(table1,['table_clustsize_totmol_totgenes_lev2_',date,'.txt']);
toc
