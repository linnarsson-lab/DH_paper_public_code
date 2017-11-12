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
strain_1 = strain_dorsalhorn_all(goodwells);

tissue_uni = unique(tissue_all_1);

% hbbgene = cellfun(@(x) ~isempty(strfind(x,'Hbb-')),geneuni(:,1)) | cellfun(@(x) ~isempty(strfind(x,'Hba-')),geneuni(:,1));
% indgood = mean(moldata_1,2)>(50/length(chip_1)) & sum(moldata_1>0,2)<0.5*length(chip_1) & ~hbbgene;
% moldata_1 = moldata_1(indgood,:);
% geneuni_1 = geneuni(indgood,:);

geneuni_1 = geneuni(:,1);
cellid_1 = cell(length(chip_1),1);
for i=1:length(chip_1)
    cellid_1{i} = [chip_1{i},'_',well_1{i}];
end


% % % % % % % % % % % % % % % % % % % % % % % %
% %Extracting neurons
% markertable = loadCellFile_turbo('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/lev1_markertable_allcells_backspinv7_lev5_20-Jun-2016.txt',1);
% markertable = markertable(2:end,1:5);
% neuronsind = find(strcmpi(markertable(:,5),'neuron'));
% cellid_neurons = markertable(neuronsind,1);
% [~,loc_neurons] = ismember(cellid_neurons, cellid_1);
% loc_neurons(loc_neurons==0) = [];
% 
% %loc_nonneurons = setdiff(cellid_1, cellid_neurons
% loc_nonneurons = setdiff([1:length(cellid_1)], loc_neurons);
% 
% 
% 
% % % % % % % % % % % % % % % %
% % t-test to exclude non-neuron genes
% [~,p] = ttest2(log2(moldata_1(:,loc_neurons)+1)',log2(moldata_1(:,loc_nonneurons)+1)','tail','right');
% %figure;plot(sort(p))
% neurons_genes_ind = intersect(setdiff([1:length(geneuni_1)],fdr_proc(1-p,0.5)),find(sum(moldata_1(:,loc_neurons),2)>10));
% %find(strcmpi(geneuni_1(neurons_genes_ind),'aqp4'))
% %find(strcmpi(geneuni_1(neurons_genes_ind),'sst'))
% % % % % % % % % % % % % % % %
% 
% 
% 
% 
% 
% 
% geneuni_1 = geneuni_1(neurons_genes_ind);
% 
% 
% moldata_1 = moldata_1(neurons_genes_ind,loc_neurons);
% chip_1 = chip_1(loc_neurons);
% well_1 = well_1(loc_neurons);
% age_1 = age_1(loc_neurons);
% sex_1 = sex_1(loc_neurons);
% diam_1 = diam_1(loc_neurons);
% tissue_all_1 = tissue_all_1(loc_neurons);
% total_mol_1 = total_mol_1(loc_neurons);
% image_all_1 = image_all_1(loc_neurons);
% moldata_spikes = moldata_spikes(:,loc_neurons);
% cellid_1 = cell(length(chip_1),1);
% for i=1:length(chip_1)
%     cellid_1{i} = [chip_1{i},'_',well_1{i}];
% end
% moldata_1 = moldata_1./repmat(total_mol_1,length(geneuni_1),1)*1e4;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% load /data2/C1_stuff/DentGyr/data_after_backspinv2_nonneurons_lev4_DG_23-Dec-2015.mat
lev2_annot = loadCellFile_turbo('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/classif_all_noIEG_bsv7_lev5_22-Jul-2016_excl.txt',1);
lev2_annot = lev2_annot(2:end,1:7);

grlev2 = lev2_annot(:,7);
% grlev2_uni = unique(grlev2);
grlev2_uni = loadCellFile('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/lev2_dendrogram_order_nov25_2016.txt');
grlev2_uni = regexprep(grlev2_uni,' ','');
grlev2_uni = flipud(grlev2_uni);
standardnames = {'Glut1','Glut2','Glut3','Glut4','Glut5','Glut6','Glut7','Glut8','Glut9','Glut10','Glut11','Glut12','Glut13','Glut14','Glut15',.....
    'Gaba1','Gaba2','Gaba3','Gaba4','Gaba5','Gaba6','Gaba7','Gaba8','Gaba9','Gaba10','Gaba11','Gaba12','Gaba13','Gaba14','Gaba15'};

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
strain_1 = strain_1(loc);


indgood = [1:length(moldata_1(1,:))];
indgood2 = indgood;
% indgood2 = find(sorted_total_mol_1>2500);
T_cells_tmp = zeros(length(grlev2),1);
grlev2_standardnames = cell(size(grlev2));
for i=1:length(grlev2_uni)
    in = strcmpi(grlev2, grlev2_uni{i});
    T_cells_tmp(in) = i;
    grlev2_standardnames(in) = repmat(standardnames(i),sum(in),1);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % %
T_cells_tmp_uni = unique(T_cells_tmp);
T_cells_tmp_uni(T_cells_tmp_uni==0) = [];

[~,xi] = sort(T_cells_tmp);
xi(T_cells_tmp(xi)==0) = [];

data_sorted_all = moldata_1(:,xi);
T_cells_tmp = T_cells_tmp(xi);
age_1 = age_1(xi);
sex_1 = sex_1(xi);
total_mol_1 = total_mol_1(xi);
cellid_1 = cellid_1(xi);
strain_1 = strain_1(xi);
grlev2 = grlev2(xi);
grlev2_standardnames = grlev2_standardnames(xi);



table1 = [cellid_1';m2c([age_1;sex_1;total_mol_1]);grlev2';grlev2_standardnames';strain_1];
table1 = [{'cellid';'age(days)';'sex(female=1)';'total mol';'cluster name';'standard_names';'strain'},table1];

saveCellFile(table1,['cell_annotation_C1_data_DHneurons_',date,'.txt']);

clear table
T = table(data_sorted_all,'RowNames',geneuni_1);
writetable(T,['moldata_C1_data_DHneurons_no_annot_',date,'.txt'],'WriteRowNames',true);


toc







