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
markertable = loadCellFile_turbo('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/lev1_markertable_allcells_backspinv7_lev5_20-Jun-2016.txt',1);
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


ind_gr_tmp_mark = [xi0(1:20,:);xi0p5(1:20,:);xi1(1:20,:)];
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


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% TSNE
data = log2_center(moldata_1(corr_filt,indgood2));

indgood = [1:length(data(1,:))];

no_dims = 2;
initial_dims = 80;length(data(:,1));%200;
perplexity = 30;
epsilon = 100;
dist_flag = 2;
max_iter = 2000;

tic
mappedX_cell = tsne((data(:,indgood))', [], no_dims, initial_dims, perplexity,epsilon,dist_flag,max_iter);
toc

       

markervec = '..........so*+><p';    
figure;
set(gcf,'color','w','position',[20,20,1000,800])   
axes('position',[0.01,0.01,0.65,0.97])
% % colors = distinguishable_colors(length(T_cells_tmp_uni));
% colors = distinguishable_colors(15);
% colors = [colors;colors];
% % colors([10,13,27,29],:) = colors([13,10,29,27],:);
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
colors = [colors_glut;color_gaba];
for i=1:length(T_cells_tmp_uni)
    in = T_cells_tmp(indgood2)==T_cells_tmp_uni(i);
    if mod(i,8)+1==1
        plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',10);hold on;
    else
        plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',10);hold on;
    end
    xcent = median(mappedX_cell(in,1));
    ycent = median(mappedX_cell(in,2));
%     htxt = text(xcent,ycent,upper(regexprep(grlev2_uni{i},'_','-')));
    htxt = text(xcent,ycent,num2str(i));
    set(htxt,'fontsize',6,'backgroundcolor',0.8*[1,1,1]);
end
axis tight
axis off
legend([upper(regexprep(grlev2_uni,'_','-'))])
title(['tsne neurons perplexity = ',num2str(perplexity),';Ngenes=',...
    num2str(length(data(:,1))),';NPCA=',num2str(initial_dims),';eps=',num2str(epsilon),';Niter=',num2str(max_iter)]);
eval(['export_fig backspinclusters_only_neurons_tsne_perplexity',num2str(perplexity),'_npca',num2str(initial_dims),'_epsilon',num2str(epsilon)...
    ,'_',date,'.pdf']);

figure;
set(gcf,'color','w','position',[20,20,1000,800])
axes('position',[0.01,0.01,0.65,0.97])
% plot(mappedX_cell(:,1),mappedX_cell(:,2),'.w'); hold on; axis tight
markervec = '...............so*+><p';
% colors = distinguishable_colors(length(T_cells_tmp_uni));
for i=1:length(T_cells_tmp_uni)
    in = T_cells_tmp==T_cells_tmp_uni(i);
    plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',1);hold on;
end
leg_str = upper(regexprep(grlev2_uni,'_','-'));
for i=1:length(T_cells_tmp_uni)
    leg_str{i} = [num2str(i),'-',leg_str{i}];
    in = T_cells_tmp==T_cells_tmp_uni(i);
    if i<=15
        text(mappedX_cell(in,1),mappedX_cell(in,2),num2str(T_cells_tmp_uni(i)),'color',colors(i,:),'fontsize',8,'fontweight','bold');hold on;
    else
        text(mappedX_cell(in,1),mappedX_cell(in,2),num2str(T_cells_tmp_uni(i)-15),'color',colors(i,:),'fontsize',8,'fontweight','bold');hold on;
    end
    xcent = median(mappedX_cell(in,1));
    ycent = median(mappedX_cell(in,2));
    %     htxt = text(xcent,ycent,upper(regexprep(grlev2_uni{i},'_','-')));
    htxt = text(xcent,ycent,num2str(i));
    set(htxt,'fontsize',6,'backgroundcolor',0.8*[1,1,1]);
end
axis tight
axis off
set(gca,'ydir','reverse');
legend(leg_str)
title(['tsne neurons perplexity = ',num2str(perplexity),';Ngenes=',...
    num2str(length(data(:,1))),';NPCA=',num2str(initial_dims),';eps=',num2str(epsilon),';Niter=',num2str(max_iter)]);
eval(['export_fig backspinclusters_withnumbers_only_neurons_tsne_perplexity',num2str(perplexity),'_npca',num2str(initial_dims),'_epsilon',num2str(epsilon)...
        ,'_',date,'.pdf']);    
    

figure;
set(gcf,'color','w','position',[20,20,600,600])
markergene = total_mol_1(indgood2);
c_rgb = [1,0,0];rand([1,3]);
markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
    interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
    ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
scatter(mappedX_cell(:,1),mappedX_cell(:,2),80,markergene_color,'.')
axis tight
title('total mol');

eval(['export_fig totalmol_only_neurons_tsne_perplexity',num2str(perplexity),'_npca',num2str(initial_dims),'_epsilon',num2str(epsilon)...
        ,'_',date,'.pdf']);


% % % % % % % % % %    
% Markers Glut1
    
    
list = {'Slc17a6','Reln','Nmur2','Tac2','Nmu','Cck','Maf','Cpne4','Trh','Meis2','Tac1','Lamp5'};
figure;
set(gcf,'color','w','position',[-1000,20,1150,800])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata_1(strcmpi(geneuni_1,genePlot),indgood2));
    if length(unique(markergene))>1
        markergene(markergene>prctile(markergene,95)) = prctile(markergene,100);
        markergene(markergene<prctile(markergene,5)) = prctile(markergene,3);
    end
    c_rgb = [1,0,0];rand([1,3]);
    %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
    %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,4,i);
    scatter(mappedX_cell(:,1),mappedX_cell(:,2),80,markergene_color,'.')
    %     set(gca,'xlim',[-150,150],'ylim',[-150,150])
    axis tight
    axis off
%     axis equal
    title(genePlot);
end


eval(['export_fig glutmarkers1_only_neurons_tsne_perplexity',num2str(perplexity),'_',date,'.pdf']);
    
% % % % % % % %
% Glut2

list = {'Slc17a6','Elavl4','Meis2','Lypd1','Nts','Qrfpr','Ecel1','Sst','Pcp4','Ucn3','Npff','Mia'};
figure;
set(gcf,'color','w','position',[20,20,1150,800])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata_1(strcmpi(geneuni_1,genePlot),indgood2));
    if length(unique(markergene))>1
        markergene(markergene>prctile(markergene,95)) = prctile(markergene,100);
        markergene(markergene<prctile(markergene,5)) = prctile(markergene,3);
    end
    c_rgb = [1,0,0];rand([1,3]);
    %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
    %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,4,i);
    scatter(mappedX_cell(:,1),mappedX_cell(:,2),80,markergene_color,'.')
    %     set(gca,'xlim',[-150,150],'ylim',[-150,150])
    axis tight
    axis off
%     axis equal
    title(genePlot);
end

eval(['export_fig glutmarkers2_only_neurons_tsne_perplexity',num2str(perplexity),'_',date,'.pdf']);



% % % % % % %
% Gaba1

list = {'Slc32a1','Gad1','Gad2','Cplx1','Crabp1','Spp1','Krt17','Sorcs3','Sst','Col25a1','Tac1','Calb2'};
figure;
set(gcf,'color','w','position',[20,20,1150,800])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata_1(strcmpi(geneuni_1,genePlot),indgood2));
    if length(unique(markergene))>1
        markergene(markergene>prctile(markergene,95)) = prctile(markergene,100);
        markergene(markergene<prctile(markergene,5)) = prctile(markergene,3);
    end
    c_rgb = [1,0,0];rand([1,3]);
    %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
    %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,4,i);
    scatter(mappedX_cell(:,1),mappedX_cell(:,2),80,markergene_color,'.')
    %     set(gca,'xlim',[-150,150],'ylim',[-150,150])
    axis tight
    axis off
%     axis equal
    title(genePlot);
end


eval(['export_fig gabamarkers1_only_neurons_tsne_perplexity',num2str(perplexity),'_',date,'.pdf']);

% % % % % %
% Gaba2   

list = {'Slc32a1','Gad1','Gad2','Gal','Pnoc','Rspo3','Npy','Qrfpr','Crabp1','Ecel1','Tac2'};
figure;
set(gcf,'color','w','position',[20,20,1150,800])
for i=1:length(list)
    genePlot = list{i};
    markergene = (moldata_1(strcmpi(geneuni_1,genePlot),indgood2));
    if length(unique(markergene))>1
        markergene(markergene>prctile(markergene,95)) = prctile(markergene,100);
        markergene(markergene<prctile(markergene,5)) = prctile(markergene,3);
    end
    c_rgb = [1,0,0];rand([1,3]);
    %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
    %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    subplot(3,4,i);
    scatter(mappedX_cell(:,1),mappedX_cell(:,2),80,markergene_color,'.')
    %     set(gca,'xlim',[-150,150],'ylim',[-150,150])
    axis tight
    axis off
%     axis equal
    title(genePlot);
end


eval(['export_fig gabamarkers2_only_neurons_tsne_perplexity',num2str(perplexity),'_',date,'.pdf']);
   
 



    
    

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% All genes in same plot - markers overlay

list = {'Tac1','Qrfpr','Pvalb','Crabp1','Spp1','Gad2','Tac2','Npy','Gal','Cplx1','Calb1','Cck'};
colors = distinguishable_colors(length(list)+1);    
[~,loc] = ismember(list,geneuni_1);
[maxval,imax] = max(2.^(log2(moldata_1(loc,indgood2)+1) - repmat(max(log2(moldata_1(loc,indgood2)+1),[],2),1,length(indgood2)) ));
imax(maxval<0.05) = length(list)+1;
markervec = '...........so*+><p';
rng(336514)
colorvec = rand(length(list)+1,3);
colorvec(end,:) = [0.8,0.8,0.8];
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=1:length(list)+1
    in = imax==i;
   if mod(i,8)+1==1
        plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',16);hold on;
    else
        plot(mappedX_cell(in,1),mappedX_cell(in,2),[markervec(mod(i,8)+1)],'color',colors(i,:),'markersize',16);hold on;
    end
end
axis tight
axis off
% axis equal
legend([list,'none'],'location','southeast') 
title(['preplexity=',num2str(perplexity),',nPCA=',num2str(initial_dims),',eps=',num2str(epsilon),',Ngenes=',num2str(length(data(:,1)))])

eval(['export_fig mixed_tsne_perplexity',num2str(perplexity),'_npca',num2str(initial_dims),'_epsilon',num2str(epsilon)...
        ,'_markersOverlay','_',date,'.pdf']);


 
% % % % % % % % % % % % % % % % % % % % % % % % % 
% T_cells_tmp = zeros(length(grlev2),1);
% for i=1:length(grlev2_uni)
%     T_cells_tmp( strcmpi(grlev2, grlev2_uni{i}) ) = i;
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% T_cells_tmp_uni = unique(T_cells_tmp);
% T_cells_tmp_uni(T_cells_tmp_uni==0) = [];
% dendrogram_order = {'GABA_Cplx1_Crabp1_Spp1pos','GABA_Cplx1_Crabp1_Spp1neg','GABA_Cplx1_Sorcs3','GABA_Cplx1_Sst',...
%     'GABA_Cplx1_Col25a1','GABA_Calb2_Cplx1','GABA_Cplx1_Adamts5','GABA_Calb2_Tac1','GABA_Calb2_Krt17','GABA_Npy_Crabp1'...
%     ,'GABA_Npy_Qrfpr','GABA_Npy_only','GABA_Tac2','GABA_Gal_Rspo3','GABA_Gal_Pnoc','Glut_Gal',...
%     'Glut_undef_GABA','Glut_Lypd1','Glut_Meis2','Glut_Elavl4_only','Glut_Qrfpr','Glut_Tac1_1','Glut_Tac1_2'...
%     ,'Glut_Reln_Npff','Glut_Reln_Nmur2','Glut_Tac2_1','Glut_Tac2_2','Glut_Tac2_3',...
%     'Glut_Nts','Glut_Cck_Trh','Glut_Cck_Cpne4','Glut_Cck_Maf'};
% [~,loc] = ismember(dendrogram_order,grlev2_uni);
% Ttmp = zeros(size(T_cells_tmp));
% for i=1:length(loc)
%     Ttmp(T_cells_tmp==loc(i)) = i;
% end
% T_cells_tmp = Ttmp;
% T_cells_tmp_uni = unique(T_cells_tmp);
% T_cells_tmp_uni(T_cells_tmp_uni==0) = [];
% [~,xi] = sort(T_cells_tmp);
% 
% data_sorted_all = moldata_1(:,xi);
% T_cells_tmp = T_cells_tmp(xi);
% data_sorted_all(:,T_cells_tmp==0) = [];
% T_cells_tmp(T_cells_tmp==0) = [];
% 
% sumpergene = sum(data_sorted_all,2);
% meanpergene = mean(data_sorted_all,2);%
% molenrich_mat = zeros(length(data_sorted_all(:,1)),length(T_cells_tmp_uni));
% meangr_mat = zeros(length(data_sorted_all(:,1)),length(T_cells_tmp_uni));
% meangrpos_mat = zeros(length(data_sorted_all(:,1)),length(T_cells_tmp_uni));
% gr_center = zeros(length(T_cells_tmp_uni),1);
% for jjj=1:length(T_cells_tmp_uni)
%     jjj
%     gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
%     molenrich_mat(:,jjj) = mean(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./meanpergene;
%     meangrpos_mat(:,jjj) = mean(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj))>0,2);
% end
% [~,grord] = sort(gr_center);
% meangrpos_mat(meangrpos_mat<0.2) = 0;
% molenrich_mat(meanpergene==0 | isnan(meanpergene),:) = 0;
% meangrpos_mat = meangrpos_mat(:,grord);
% molenrich_mat = molenrich_mat(:,grord);
% [~,xi0] = sort(molenrich_mat.*meangrpos_mat.^0.001,'descend');
% [~,xi0p5] = sort(molenrich_mat.*meangrpos_mat.^0.5,'descend');
% [~,xi1] = sort(molenrich_mat.*meangrpos_mat.^1,'descend');
% 
% 
% % ind_gr_tmp_mark = [xi0(1:3,:);xi0p5(1:3,:);xi1(1:3,:)];
% ind_gr_tmp_mark = [xi0p5(1:3,:)];
% ind_gr_tmp_mark = flipud(ind_gr_tmp_mark(:));
% rmv = false(size(ind_gr_tmp_mark));
% for i=2:length(ind_gr_tmp_mark)
%     if ismember(ind_gr_tmp_mark(i),ind_gr_tmp_mark(1:i-1))
%         rmv(i) = true;
%     end
% end
% ind_gr_tmp_mark(rmv) = [];
% gr_tmp_mark = geneuni_1(ind_gr_tmp_mark);
% gr_tmp_mark = (gr_tmp_mark(:));
% molenrich_mat_mark = zeros(length(ind_gr_tmp_mark),length(T_cells_tmp_uni));
% meangrpos_mat_mark = zeros(length(ind_gr_tmp_mark),length(T_cells_tmp_uni));
% order_gr = [T_cells_tmp(diff(T_cells_tmp)~=0)',T_cells_tmp(end)];
% for jjj=1:length(T_cells_tmp_uni)
%     jjj
%     molenrich_mat_mark(:,jjj) = mean(data_sorted_all(ind_gr_tmp_mark,T_cells_tmp==order_gr(jjj)),2)./meanpergene(ind_gr_tmp_mark);
%     meangrpos_mat_mark(:,jjj) = mean(data_sorted_all(ind_gr_tmp_mark,T_cells_tmp==order_gr(jjj))>0,2);
% end
% [~,imax] = max(molenrich_mat_mark.*meangrpos_mat_mark.^0.5,[],2);
% [~,xi] = sort(imax);
% ind_gr_tmp_mark = [ind_gr_tmp_mark(xi)];
% 
% datamarkers = data_sorted_all(ind_gr_tmp_mark,:);
% datamarkers_cn = cent_norm(log2(datamarkers+1));
% gr_tmp_mark = [gr_tmp_mark(xi)];
% 
% figure;
% set(gcf,'position',[100,100,700,900],'color','w')
% ax1 = axes('position',[0.1,0.1,0.88,0.76]);
% imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),3),prctile(datamarkers_cn(:),97)]);
% hold on;
% linewid =0.5;
% bor_color = 'grey11';%'green1';%
% cells_bor = find(diff(T_cells_tmp)~=0);
% for jj=1:length(cells_bor)
%     plot(cells_bor(jj)*[1,1]-0.5,[1,length(gr_tmp_mark)],'-','linewidth',linewid,'color',get_RGB(bor_color))
% end
% set(gca,'xtick',gr_center,'xticklabel',regexprep(dendrogram_order,'_','-'),'XTickLabelRotation',45,'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 8)
% colormap('summer');
% freezeColors(gca);
% axes('position',[0.1,0.87,0.88,0.01]);
% % bar(data_sorted_all(strcmpi(geneuni_1,'slc32a1'),:)); axis tight
% hb1 = bar(data_sorted_all(strcmpi(geneuni_1,'gad1'),:)>0); axis tight;
% set(hb1,'facecolor','r');
% set(gca,'xtick',[],'ytick',[0.5],'yticklabel','Gad1')
%    
% axes('position',[0.1,0.88,0.88,0.01]);
% % bar(data_sorted_all(strcmpi(geneuni_1,'slc32a1'),:)); axis tight
% hb2 = bar(data_sorted_all(strcmpi(geneuni_1,'slc17a6'),:)>0); axis tight
% set(hb2,'facecolor','b');
% set(gca,'xtick',[],'ytick',[0.5],'yticklabel','Slc17a6')
% 
% 
% eval(['export_fig DH_neurons_86x1621_markertable_dendrogramOrder','_',date,'.pdf']);
%     
%     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    