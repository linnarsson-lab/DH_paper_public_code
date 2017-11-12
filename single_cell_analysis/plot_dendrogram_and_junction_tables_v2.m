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

hbbgene = cellfun(@(x) ~isempty(strfind(x,'Hbb-')),geneuni(:,1)) | cellfun(@(x) ~isempty(strfind(x,'Hba-')),geneuni(:,1));
indgood = mean(moldata_1,2)>(50/length(chip_1)) & sum(moldata_1>0,2)<0.5*length(chip_1) & ~hbbgene;
moldata_1 = moldata_1(indgood,:);
geneuni_1 = geneuni(indgood,:);
cellid_1 = cell(length(chip_1),1);
for i=1:length(chip_1)
    cellid_1{i} = [chip_1{i},'_',well_1{i}];
end


% % % % % % % % % % % % % % % % % % % % % % %
%Extracting neurons
markertable = loadCellFile_turbo('/data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/lev1_markertable_allcells_backspinv7_lev5_20-Jun-2016.txt',1);
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
neurons_genes_ind = intersect(setdiff([1:length(geneuni_1)],fdr_proc(1-p,0.5)),find(sum(moldata_1(:,loc_neurons),2)>10));
%find(strcmpi(geneuni_1(neurons_genes_ind),'aqp4'))
%find(strcmpi(geneuni_1(neurons_genes_ind),'sst'))
% % % % % % % % % % % % % % %






geneuni_1 = geneuni_1(neurons_genes_ind);


moldata_1 = moldata_1(neurons_genes_ind,loc_neurons);
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
lev2_annot = loadCellFile_turbo('/data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/classif_all_noIEG_bsv7_lev5_22-Jul-2016_excl.txt',1);
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


indgood = [1:length(moldata_1(1,:))];
indgood2 = indgood;
% indgood2 = find(sorted_total_mol_1>2500);
T_cells_tmp = zeros(length(grlev2),1);
for i=1:length(grlev2_uni)
    T_cells_tmp( strcmpi(grlev2, grlev2_uni{i}) ) = i;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % %
T_cells_tmp_uni = unique(T_cells_tmp);
T_cells_tmp_uni(T_cells_tmp_uni==0) = [];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

log_m_gr_lev2 = zeros(length(geneuni_1),length(T_cells_tmp_uni));
se_log_m_gr_lev2 = zeros(length(geneuni_1),length(T_cells_tmp_uni));
log_m_gr_lev2_norm = zeros(length(geneuni_1),length(T_cells_tmp_uni));
se_log_m_gr_lev2_norm = zeros(length(geneuni_1),length(T_cells_tmp_uni));
present_mean = zeros(length(geneuni_1),length(T_cells_tmp_uni));
moldata_1_norm = moldata_1./repmat(sum(moldata_1),length(geneuni_1),1)*1e4;
for i=1:length(T_cells_tmp_uni)
    i
    in = T_cells_tmp==T_cells_tmp_uni(i);
    log_m_gr_lev2(:,i) = mean(log2(moldata_1(:,in)+1),2);
    se_log_m_gr_lev2(:,i) = std(log2(moldata_1(:,in)+1),[],2)/sqrt(sum(in));
    log_m_gr_lev2_norm(:,i) = mean(log2(moldata_1_norm(:,in)+1),2);
    se_log_m_gr_lev2_norm(:,i) = std(log2(moldata_1_norm(:,in)+1),[],2)/sqrt(sum(in));
    present_mean(:,i) = mean(moldata_1(:,in)>0,2);
end


mean_exp_gr_lev2 = zeros(length(geneuni_1),length(T_cells_tmp_uni));
mean_exp_gr_lev2_norm = zeros(length(geneuni_1),length(T_cells_tmp_uni));
for i=1:length(T_cells_tmp_uni)
    in = T_cells_tmp==T_cells_tmp_uni(i);
    mean_exp_gr_lev2(:,i) = mean(moldata_1(:,in),2);
    mean_exp_gr_lev2_norm(:,i) = mean(moldata_1_norm(:,in),2);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % %

% sumpergene = sum(dataout_sorted_all_nonn,2);
meanpergene = mean(moldata_1,2);%
molenrich_mat = zeros(length(moldata_1(:,1)),length(T_cells_tmp_uni));
% meangr_mat = zeros(length(moldata_1(:,1)),length(T_cells_tmp_uni));
meangrpos_mat = zeros(length(moldata_1(:,1)),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
    molenrich_mat(:,jjj) = mean(moldata_1(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./meanpergene;
    meangrpos_mat(:,jjj) = mean(moldata_1(:,T_cells_tmp==T_cells_tmp_uni(jjj))>0,2);
end
[~,grord] = sort(gr_center);
meangrpos_mat = meangrpos_mat(:,grord);
molenrich_mat = molenrich_mat(:,grord);
[~,xi0] = sort(molenrich_mat.*meangrpos_mat.^0.00001,'descend');
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
corr_filt = ind_gr_tmp_mark;


% % % % % % % % % % % % % % % % % % % % % % % % % 

log_m_mat_link = log_m_gr_lev2_norm(ind_gr_tmp_mark,:);%log2(mean_exp_gr_lev2_norm(ind_gr_tmp_mark,:)+1);%
gene_link = gr_tmp_mark;
Z = linkage(log_m_mat_link','average','correlation');
Rcelltypes = corr_mat(log_m_mat_link);
D = pdist(log_m_mat_link','correlation');
leaforder = optimalleaforder(Z,D);
% {'PPR','OPC', 'COP', 'NFOL1','NFOL2','MFOL1','MFOL2','MOL1','MOL2','MOL4','MOL3','MOL5','MOL6'}
% [~,leaforder] = ismember({'PPR','OPC', 'COP', 'NFOL1','NFOL2','MFOL1','MFOL2','MOL1','MOL4','MOL2','MOL3','MOL5','MOL6'},lev2_gr_uni);
Z = linkage(log_m_mat_link(:,leaforder)','average','correlation');


figure('color','w','position',[200,200,1450,800],'visible','on');
axes('position',[0.03,0.05,0.4,0.9])
[H,T,outperm] = dendrogram(Z,length(T_cells_tmp_uni),'Orientation','left','reorder',...
    [1:length(T_cells_tmp_uni)],'labels',regexprep(grlev2_uni(leaforder),'_','-'),'colorthreshold',0.2);
set(H,'linewidth',2)
xl = get(gca,'xlim');
set(gca,'ylim',[0.5,length(T_cells_tmp_uni)+0.5],'xtick',[],'ydir','normal','fontsize',8,'xlim',[0,0.7])
% set(gca,'xlim',[0.5,length(T_cells_tmp_uni)+1],'ytick',[],'fontsize',8)
% rotateticklabel(gca,45);

axes('position',[0.58,0.05,0.35,0.9])
imagesc(Rcelltypes(leaforder,leaforder),[0.3,1]) 
set(gca,'xtick',[],'xticklabel',leaforder,'xdir','normal','ydir','normal','ytick',[])
colormap('summer')
hcm = colorbar;
set(hcm,'position',[0.96,0.7,0.02,0.25],'ytick',[-0.4:0.1:1],'ydir','normal')
% eval(['export_fig oligos_dendogram_by_lev2_',date,'.pdf'])
save2pdf(['DH_dendogram_by_lev2_',date,'.pdf'],gcf,300)

figure('color','w','position',[200,200,800,1000],'visible','on');
axes('position',[0.03,0.03,0.85,0.96])
[H,T,outperm] = dendrogram(Z,length(T_cells_tmp_uni),'Orientation','left','reorder',...
    [1:length(T_cells_tmp_uni)],'labels',regexprep(grlev2_uni(leaforder),'_','-'),'colorthreshold',0.2);
set(H,'linewidth',1)
for i=1:length(H)
    xd = get(H(i),'XData');
    yd = get(H(i),'YData');
    text(xd(2),mean(yd(2:3)),num2str(i),'color','r');
end
eval(['export_fig -r 2000 DH_dendogram_by_lev2_withnumbers_',date,'.pdf'])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


n = length(grlev2_uni);
clut2leaf = cell(max(Z(1:2*length(Z))),1);
clut2leaf(1:n) = m2c(1:n)';
clut2branch = cell(max(Z(1:2*length(Z)))+1,1);
clutaveR(1:length(Rcelltypes)) = 1;
clutham(1:length(Rcelltypes)) = 1;
bin_pat_3 = zeros(size(log_m_mat_link));
for i=1:length(Z)
    clut2leaf{i+n} = [clut2leaf{Z(i,1)};clut2leaf{Z(i,2)}];
    clut2branch{Z(i,1)} = i;
    clut2branch{Z(i,2)} = i;
    bin_pat_3(:,i+length(Rcelltypes)) = mean(bin_pat_3(:,Z(i,1:2)),2);
    clutaveR(i+length(Rcelltypes)-1) = corr(bin_pat_3(:,Z(i,1)),bin_pat_3(:,Z(i,2)));
    clutham(i+length(Rcelltypes)-1) = sum(bin_pat_3(:,Z(i,1))==bin_pat_3(:,Z(i,2)))/length(bin_pat_3(:,1));
end

binary_profiles = zeros(n,length(clut2leaf));
binary_profiles_branch = zeros(1,length(clut2leaf));
for i=1:length(binary_profiles(1,:))-1
    binary_profiles(clut2leaf{i},i) = 1;
    binary_profiles_branch(i) = clut2branch{i};
end
binary_profiles_names = cell(size(binary_profiles));
for i=1:length(binary_profiles_names)
    binary_profiles_names(1:sum(binary_profiles(:,i)>0),i) = grlev2_uni(leaforder(binary_profiles(:,i)>0));
end

binary_profiles_sc = false(length(geneuni_1), length(binary_profiles_names(1,:)));
for i=1:length(binary_profiles_names(1,:))
    in = find(~cellfun(@isempty, binary_profiles_names(:,i)));
    for j=1:length(in)
        binary_profiles_sc(strcmpi(grlev2,binary_profiles_names{in(j),i}),i) = true;
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % 


ind_branch_sep = zeros(length(geneuni_1),max(binary_profiles_branch));%cell(100,1);
pcomp1_all = ind_branch_sep;
qcomp1_all = ind_branch_sep;
pcomp2_all = ind_branch_sep;
qcomp2_all = ind_branch_sep;


frac_gr1_mat = zeros(length(geneuni_1),max(binary_profiles_branch));
frac_gr2_mat = zeros(length(geneuni_1),max(binary_profiles_branch));
names_cell = cell(max(binary_profiles_branch),2);
for i=1:max(binary_profiles_branch)
    i
    compind = find(binary_profiles_branch==i);
    gr1 = binary_profiles_sc(:,compind(1));
    gr2 = binary_profiles_sc(:,compind(2));
    pre_frac_gr1 = mean(moldata_1(:,gr1)>0,2);
    pre_frac_gr2 = mean(moldata_1(:,gr2)>0,2);
    pre_frac_all = mean(moldata_1>0,2);
    gr1_name = binary_profiles_names(:,compind(1));
    gr1_num = find(binary_profiles(:,compind(1)));
    gr2_name = binary_profiles_names(:,compind(2));
    gr2_num = find(binary_profiles(:,compind(2)));
    gr1_name(cellfun(@isempty,gr1_name)) = [];
    gr2_name(cellfun(@isempty,gr2_name)) = [];
    tmp = [];
    tmpnum = [];
    for j=1:length(gr1_name)
        tmp = [tmp,',',gr1_name{j}];
        tmpnum = [tmpnum,',',num2str(gr1_num(j))];
    end
    gr1_name = ['(',tmpnum,')',',',tmp];
    tmp = [];
    tmpnum = [];
    for j=1:length(gr2_name)
        tmp = [tmp,',',gr2_name{j}];
        tmpnum = [tmpnum,',',num2str(gr2_num(j))];
    end
    gr2_name = tmp;
    gr2_name = ['(',tmpnum,')',',',tmp];
    names_cell(i,:) = [{gr1_name},{gr2_name}];
    score_comp1 = pre_frac_gr1-pre_frac_gr2;
    pbino_comp1 = binocdf(sum(gr1)*pre_frac_gr1,sum(gr1),pre_frac_all,'upper');
    qbino_comp1 = qval_from_pval(pbino_comp1);
    score_comp2 = pre_frac_gr2-pre_frac_gr1;
    pbino_comp2 = binocdf(sum(gr2)*pre_frac_gr2,sum(gr2),pre_frac_all,'upper');
    qbino_comp2 = qval_from_pval(pbino_comp2);
    
    [~,xi1] = sort(score_comp1,'descend');
    [~,xi2] = sort(score_comp2,'descend');
    ind_branch_sep(xi1(1:50),i) = score_comp1(xi1(1:50));
    ind_branch_sep(xi2(1:50),i) = -score_comp2(xi2(1:50));
    pcomp1_all(:,i) = pbino_comp1;
%     qcomp1_all(xi1(1:50),i) = qbino_comp1(xi1(1:50));
    pcomp2_all(:,i) = pbino_comp2;
%     qcomp2_all(xi2(1:50),i) = qbino_comp2(xi2(1:50));
        
    frac_gr1_mat(:,i) = pre_frac_gr1;
    frac_gr2_mat(:,i) = pre_frac_gr2;
end
for i=1:length(pcomp1_all(:,1))
    qcomp1_all(i,:) = qval_from_pval(pcomp1_all(i,:));
    qcomp2_all(i,:) = qval_from_pval(pcomp2_all(i,:));
end


table1 = cell(53*max(binary_profiles_branch),12);
k=0;
for i=1:max(binary_profiles_branch); 
    i, 
    in1 = find(ind_branch_sep(:,i)>0);
    in2 = find(ind_branch_sep(:,i)<0);
    [~,xi1] = sort(ind_branch_sep(in1,i),'descend');
    [~,xi2] = sort(-ind_branch_sep(in2,i),'descend');
    in1 = in1(xi1);
    in2 = in2(xi2);
    tmp = [geneuni_1(in1), m2c([ind_branch_sep(in1,i),pcomp1_all(in1,i),qcomp1_all(in1,i),frac_gr1_mat(in1,i),frac_gr2_mat(in1,i)]),...
        geneuni_1(in2), m2c([ind_branch_sep(in2,i),pcomp2_all(in2,i),qcomp2_all(in2,i),frac_gr1_mat(in2,i),frac_gr2_mat(in2,i)])];
%     tmp = [ [{['junction #',num2str(i)]},cell(1,7)];{'left marker','left score','left mean','right mean','right marker','right score','left mean','right mean'};tmp];
    tmp = [ [{['junction #',num2str(i),',left=',names_cell{i,1},',right=',names_cell{i,2}]},cell(1,11)];...
        {'left marker','left score','left p(binomial)','left q(binomial)','left mean','right mean','right marker'...
        ,'right score','right p(binomial)','right q(binomial)','left mean','right mean'};tmp];
    table1(k+1:k+52,:) = tmp;
    k = k+53;
end

saveCellFile(table1,['DH_neurons_dendrogram_junction_split_markers_',date,'.txt']);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

ind_branch_sep2 = zeros(length(geneuni_1),max(binary_profiles_branch));%cell(100,1);
pRankSum1_all = ind_branch_sep2;
pRankSum2_all = ind_branch_sep2;
m_gr1_mat = zeros(length(geneuni_1),max(binary_profiles_branch));
m_gr2_mat = zeros(length(geneuni_1),max(binary_profiles_branch));
for i=1:max(binary_profiles_branch)
    i
    compind = find(binary_profiles_branch==i);
    gr1 = binary_profiles_sc(:,compind(1));
    gr2 = binary_profiles_sc(:,compind(2));
%     pre_frac_gr1 = mean(moldata_1(:,gr1)>0,2);
%     pre_frac_gr2 = mean(moldata_1(:,gr2)>0,2);
%     pre_frac_all = mean([pre_frac_gr1,pre_frac_gr2],2);
%     pre_frac_gr1 = pre_frac_gr1./pre_frac_all;
%     pre_frac_gr2 = pre_frac_gr2./pre_frac_all;
    gr1_name = binary_profiles_names(:,compind(1));
    gr1_num = find(binary_profiles(:,compind(1)));
    gr2_name = binary_profiles_names(:,compind(2));
    gr2_num = find(binary_profiles(:,compind(2)));
    gr1_name(cellfun(@isempty,gr1_name)) = [];
    gr2_name(cellfun(@isempty,gr2_name)) = [];
    tmp = [];
    tmpnum = [];
    for j=1:length(gr1_name)
        tmp = [tmp,',',gr1_name{j}];
        tmpnum = [tmpnum,',',num2str(gr1_num(j))];
    end
    gr1_name = ['(',tmpnum,')',',',tmp];
    tmp = [];
    tmpnum = [];
    for j=1:length(gr2_name)
        tmp = [tmp,',',gr2_name{j}];
        tmpnum = [tmpnum,',',num2str(gr2_num(j))];
    end
    gr2_name = tmp;
    gr2_name = ['(',tmpnum,')',',',tmp];
    names_cell(i,:) = [{gr1_name},{gr2_name}];
    m_gr1 = mean(log2(moldata_1(:,gr1)+1),2);
    m_gr2 = mean(log2(moldata_1(:,gr2)+1),2);
    m_all = mean([m_gr1,m_gr2],2);
    p_rs1 = nan(length(geneuni_1),1);
    p_rs2 = nan(length(geneuni_1),1);
    for j=1:length(geneuni_1);
        p_rs1(j) = ranksum(moldata_1(j,gr1),moldata_1(j,gr2),'tail','right');
        p_rs2(j) = ranksum(moldata_1(j,gr1),moldata_1(j,gr2),'tail','left');
    end
    q1 = qval_from_pval(p_rs1);
    q2 = qval_from_pval(p_rs2);
    
    score_comp1 = m_gr1-m_gr2;
    score_comp2 = m_gr2-m_gr1;    
    [~,xi1] = sort(score_comp1,'descend');
    [~,xi2] = sort(score_comp2,'descend');
    ind_branch_sep2(xi1(1:50),i) = score_comp1(xi1(1:50));
    ind_branch_sep2(xi2(1:50),i) = -score_comp2(xi2(1:50));
    pRankSum1_all(:,i) = p_rs1;
    pRankSum2_all(:,i) = p_rs2;
    
    m_gr1_mat(:,i) = m_gr1;
    m_gr2_mat(:,i) = m_gr2;
end
for i=1:length(pRankSum1_all(:,1))
    qRankSum1_all(i,:) = qval_from_pval(pRankSum1_all(i,:));
    qRankSum2_all(i,:) = qval_from_pval(pRankSum2_all(i,:));
end


table1 = cell(53*max(binary_profiles_branch),12);
k=0;
for i=1:max(binary_profiles_branch); 
    i, 
    in1 = find(ind_branch_sep2(:,i)>0);
    in2 = find(ind_branch_sep2(:,i)<0);
    [~,xi1] = sort(ind_branch_sep2(in1,i),'descend');
    [~,xi2] = sort(-ind_branch_sep2(in2,i),'descend');
    in1 = in1(xi1);
    in2 = in2(xi2);
    tmp = [geneuni_1(in1), m2c([ind_branch_sep2(in1,i),pRankSum1_all(in1,i),qRankSum1_all(in1,i),m_gr1_mat(in1,i),m_gr2_mat(in1,i)]),...
        geneuni_1(in2), m2c([ind_branch_sep2(in2,i),pRankSum2_all(in2,i),qRankSum2_all(in2,i),m_gr1_mat(in2,i),m_gr2_mat(in2,i)])];
%     tmp = [ [{['junction #',num2str(i)]},cell(1,7)];{'left marker','left score','left mean','right mean','right marker','right score','left mean','right mean'};tmp];
    tmp = [ [{['junction #',num2str(i),',left=',names_cell{i,1},',right=',names_cell{i,2}]},cell(1,11)];...
        {'left marker','left score','left p(ranksum)','left q(ranksum)','left mean','right mean'...
        ,'right marker','right score','left p(ranksum)','left q(ranksum)','left mean','right mean'};tmp];
    table1(k+1:k+52,:) = tmp;
    k = k+53;
end

saveCellFile(table1,['DH_neurons_dendrogram_junction_split_markers_by_average_',date,'.txt']);



% % % % % % % % % % % % % % % % % % % % %
% sumpergene = sum(dataout_sorted_all_nonn,2);
meanpergene = mean(moldata_1,2);%
molenrich_mat = zeros(length(moldata_1(:,1)),length(T_cells_tmp_uni));
% meangr_mat = zeros(length(moldata_1(:,1)),length(T_cells_tmp_uni));
meangrpos_mat = zeros(length(moldata_1(:,1)),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
    molenrich_mat(:,jjj) = mean(moldata_1(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./meanpergene;
    meangrpos_mat(:,jjj) = mean(moldata_1(:,T_cells_tmp==T_cells_tmp_uni(jjj))>0,2);
    gr_size(jjj) = sum(T_cells_tmp==T_cells_tmp_uni(jjj));
end
[~,grord] = sort(gr_center);
% meangrpos_mat = meangrpos_mat(:,grord);
% molenrich_mat = molenrich_mat(:,grord);
score_mat_p1 = (meangrpos_mat.^1).*molenrich_mat;
score_mat_p1 = score_mat_p1(:,leaforder);

table1 = [ [{'clustername';'clustersize'};geneuni_1], [grlev2_uni(leaforder)'; m2c(gr_size(leaforder));  m2c(score_mat_p1)] ];
saveCellFile(table1,['DH_neurons_enrichment_score_lev2_',date,'.txt'])





% % % % % % % % % % % % % % % % % % % % % % % % % save data for cef file 
clear table
T = table([moldata_1],'RowNames',geneuni_1);
writetable(T,['moldate_DHneurons_lev2_no_annot_',date,'.txt'],'WriteRowNames',true);

annottable = [cellid_1';grlev2';m2c(sum(moldata_1))];
saveCellFile(annottable,['annotation_moldate_DHneurons_lev2_',date,'.txt'])


