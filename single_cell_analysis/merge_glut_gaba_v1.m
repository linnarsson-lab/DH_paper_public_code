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
% 
% hbbgene = cellfun(@(x) ~isempty(strfind(x,'Hbb-')),geneuni(:,1)) | cellfun(@(x) ~isempty(strfind(x,'Hba-')),geneuni(:,1));
% indgood = mean(moldata_1,2)>(50/length(chip_1)) & sum(moldata_1>0,2)<0.5*length(chip_1) & ~hbbgene;
% moldata_1 = moldata_1(indgood,:);
geneuni_1 = geneuni;

list_markers = loadCellFile('/data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/woIEGmarkers.txt');
[~,ind_markers] = ismember(list_markers,geneuni_1);

load /data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/data_glut_after_backspinv7_noIEG_lev5_dorsalhorn_glutonly_26-Jul-2016.mat
cells_sorted_glut = cells_sorted;
cells_bor_glut = cells_bor_level{end};
% 
% sorted_strain_glut = strain_1(loc);
% sorted_strain_glut_uni = unique(sorted_strain_glut);
% sorted_strain_num_glut = zeros(size(chip_1));
% for i=1:length(sorted_strain_glut_uni)
%     sorted_strain_num_glut(strcmpi(sorted_strain_glut,sorted_strain_glut_uni{i})) = i;
% end

load /data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/data_gaba_noIEG_after_backspinv7_lev5_dorsalhorn_gabaonly_26-Jul-2016.mat
cells_sorted_gaba = cells_sorted;
cells_bor_gaba = cells_bor_level{end};
[~,loc_glut_gaba] = ismember([cells_sorted_glut,cells_sorted_gaba], cellid_1);
cells_bor = [cells_bor_glut; (cells_bor_gaba+length(cells_sorted_glut))];


moldata_markers = moldata_1(:,loc_glut_gaba);

% % % % % % % % % % % % % % % % % % % % % % % % % % 
ind_gr_tmp_mark = ind_markers;


figure;
set(gcf,'position',[100,100,450,750],'color','w')

barhight = 0.89/length(ind_gr_tmp_mark);
for jj=1:length(ind_gr_tmp_mark)
    disp([num2str(jj),','])
    maxval = prctile(moldata_markers(ind_gr_tmp_mark(jj),:),99.99);%0.9*max(data_tmpgr(ind_gr_tmp_mark(jj),:));
    n = length(loc_glut_gaba);
    axes('position',[0.12,0.01+barhight*(jj-1),0.8,barhight])
    %     if mod(jj,2)==1;
    % bar(n/2,maxval,'facecolor',get_RGB('grey81'),'edgecolor','none','BarWidth',n);hold on;axis tight;axis off;
    %     end
    if mod(jj,2)==1;
        h = bar(moldata_markers(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',1);hold on;axis tight;axis off;
    else
        h = bar(moldata_markers(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',1,'facecolor',[200,0,0]/256);hold on;axis tight;axis off;
    end
    for jjj=1:length(cells_bor)
        plot(cells_bor(jjj)*[1,1]-1,[0,maxval],'--','linewidth',0.5,'color',[200,200,200]/256); hold on;
    end
    
    set(gca,'xlim',[0,n],'ylim',[0,maxval], 'fontsize', 6)
    text(-0.5,maxval/2,list_markers{jj},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',7)
end


% 
% axes('position',[0.12,0.98,0.88,0.01])
% imagesc(cent_norm(log2(sorted_strain_num)))
% set(gca,'xtick',[],'ytick',1,'yticklabel','strain')
% colormap('colorcube');freezeColors(gca);


% 
% eval(['export_fig barplots_allmarkers_noIEG_bsv7_lev',num2str(splitlev),'_',date,'.pdf'])
% 



% % % % % % % % % % % % % %
% Neurotransmitters & Receptors

list_markers = loadCellFile('/data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/neurotransm_receptors.txt');
[~,ind_markers] = ismember(list_markers,geneuni_1);

list_markers(ind_markers==0)=[];
ind_markers(ind_markers==0)=[];

ind_gr_tmp_mark = ind_markers;

figure;
set(gcf,'position',[100,100,450,750],'color','w')

barhight = 0.89/length(ind_gr_tmp_mark);
for jj=1:length(ind_gr_tmp_mark)
    disp([num2str(jj),','])
    maxval = prctile(moldata_markers(ind_gr_tmp_mark(jj),:),99.99);%0.9*max(data_tmpgr(ind_gr_tmp_mark(jj),:));
    n = length(loc_glut_gaba);
    axes('position',[0.12,0.01+barhight*(jj-1),0.8,barhight])
    %     if mod(jj,2)==1;
    % bar(n/2,maxval,'facecolor',get_RGB('grey81'),'edgecolor','none','BarWidth',n);hold on;axis tight;axis off;
    %     end
    if mod(jj,2)==1;
        h = bar(moldata_markers(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',1);hold on;axis tight;axis off;
    else
        h = bar(moldata_markers(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',1,'facecolor',[200,0,0]/256);hold on;axis tight;axis off;
    end
    for jjj=1:length(cells_bor)
        plot(cells_bor(jjj)*[1,1]-1,[0,maxval],'--','linewidth',0.5,'color',[200,200,200]/256); hold on;
    end
    
    set(gca,'xlim',[0,n],'ylim',[0,maxval], 'fontsize', 5)
    text(-0.5,maxval/2,list_markers{jj},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',6)
end


% 
% eval(['export_fig barplots_neurotransm_noIEG_bsv7_lev',num2str(splitlev),'_',date,'.pdf'])



% % % % % % % % % % % % % %
% Peptides & Receptors

list_markers = loadCellFile('/data2/C1_stuff/Dorsal_horn_MH/Analysis/analysis_jun22/peptides_receptors.txt');
[~,ind_markers] = ismember(list_markers,geneuni_1);

list_markers(ind_markers==0)=[];
ind_markers(ind_markers==0)=[];

ind_gr_tmp_mark = ind_markers;

figure;
set(gcf,'position',[100,100,450,750],'color','w')

barhight = 0.89/length(ind_gr_tmp_mark);
for jj=1:length(ind_gr_tmp_mark)
    disp([num2str(jj),','])
    maxval = prctile(moldata_markers(ind_gr_tmp_mark(jj),:),99.99);%0.9*max(data_tmpgr(ind_gr_tmp_mark(jj),:));
    n = length(loc_glut_gaba);
    axes('position',[0.12,0.01+barhight*(jj-1),0.8,barhight])
    %     if mod(jj,2)==1;
    % bar(n/2,maxval,'facecolor',get_RGB('grey81'),'edgecolor','none','BarWidth',n);hold on;axis tight;axis off;
    %     end
    if mod(jj,2)==1;
        h = bar(moldata_markers(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',1);hold on;axis tight;axis off;
    else
        h = bar(moldata_markers(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',1,'facecolor',[200,0,0]/256);hold on;axis tight;axis off;
    end
    for jjj=1:length(cells_bor)
        plot(cells_bor(jjj)*[1,1]-1,[0,maxval],'--','linewidth',0.5,'color',[200,200,200]/256); hold on;
    end
    
    set(gca,'xlim',[0,n],'ylim',[0,maxval], 'fontsize', 6)
    text(-0.5,maxval/2,list_markers{jj},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',7)
end



% 
% eval(['export_fig barplots_peptides_noIEG_bsv7_lev',num2str(splitlev),'_',date,'.pdf'])





% %% % % % % % %
% markertable



table1 = m2c([T_cells_tmp;datamarkers]');
table1 = [[{'','tissue','cluster'},gr_tmp_mark'];[cells_sorted',sorted_tissue_1',table1]];


saveCellFile(table1,['markertable_gaba_noIEG_backspinv7_lev',num2str(splitlev),'_',date,'.txt'])



