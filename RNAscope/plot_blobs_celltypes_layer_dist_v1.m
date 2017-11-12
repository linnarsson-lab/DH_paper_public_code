tic
clear all
close all

load /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/aggregate_Glut_RNAscope_counts.mat
%dataglut_LSC dataglut_CSC imagesource_glut_LSC imagesource_glut_CSC
%table_header_Glut
load /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/aggregate_GABA_RNAscope_counts.mat
%datagaba_LSC datagaba_CSC imagesource_gaba_LSC imagesource_gaba_CSC
%table_header_GABA

combi_needed = loadCellFile('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/celltypes_combinations_to_output_v2.txt');
rmv = strcmpi(combi_needed(:,3),'-') | ~cellfun(@isempty, strfind(combi_needed(:,3),'alt_'));
combi_needed(rmv,:) = [];

load ref_LSC;
load ref_CSC;
load layers_mask_LSC_ref
load layers_mask_CSC_ref

layercounts_LSC = zeros(length(combi_needed(:,1)),8);
layercounts_CSC = zeros(length(combi_needed(:,1)),8);
th = 3;
for i=1:length(combi_needed(:,1))
    i
    isgaba = ~isempty(strfind(combi_needed{i,1},'GAB'));
    isglut = ~isempty(strfind(combi_needed{i,1},'Glu'));
    tmp = strsplit(combi_needed{i,2},',');
    tmpfname = combi_needed{i,2};
    tmpfname = regexprep(tmpfname,'+','pos_');
    tmpfname = regexprep(tmpfname,'-','neg_');
    tmpfname = regexprep(tmpfname,',','-');
    
    if isgaba
        in_LSC = true(length(datagaba_LSC(:,1)),1);
        in_CSC = true(length(datagaba_CSC(:,1)),1);
        in_LSC(cellfun(@isempty, strfind(imagesource_gaba_LSC, combi_needed{i,1}))) = false;
        in_CSC(cellfun(@isempty, strfind(imagesource_gaba_CSC, combi_needed{i,1}))) = false;
        for j=1:length(tmp)
            if tmp{j}(1)=='+'
                in_LSC = in_LSC & datagaba_LSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th;
                in_CSC = in_CSC & datagaba_CSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th;
            elseif tmp{j}(1)=='-'
                in_LSC = in_LSC & ~(datagaba_LSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th);
                in_CSC = in_CSC & ~(datagaba_CSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th);
            end
        end
        
        [N,edges] = histcounts(datagaba_LSC(in_LSC,1),[1:9]);
        layercounts_LSC(i,:) = N;
        [N,edges] = histcounts(datagaba_CSC(in_CSC,1),[1:9]);
        layercounts_CSC(i,:) = N;
        
        
    elseif isglut
        in_LSC = true(length(dataglut_LSC(:,1)),1);
        in_CSC = true(length(dataglut_CSC(:,1)),1);
        in_LSC(cellfun(@isempty, strfind(imagesource_glut_LSC, combi_needed{i,1}))) = false;
        in_CSC(cellfun(@isempty, strfind(imagesource_glut_CSC, combi_needed{i,1}))) = false;
        for j=1:length(tmp)
            if tmp{j}(1)=='+'
                in_LSC = in_LSC & dataglut_LSC(:,find(strcmpi(table_header_Glut,tmp{j}(2:end)) ) )>th;
                in_CSC = in_CSC & dataglut_CSC(:,find(strcmpi(table_header_Glut,tmp{j}(2:end)) ) )>th;
            elseif tmp{j}(1)=='-'
                in_LSC = in_LSC & ~(dataglut_LSC(:,find(strcmpi(table_header_Glut,tmp{j}(2:end)) ) )>th);
                in_CSC = in_CSC & ~(dataglut_CSC(:,find(strcmpi(table_header_Glut,tmp{j}(2:end)) ) )>th);
            end
        end
        [N,edges] = histcounts(dataglut_LSC(in_LSC,1),[1:9]);
        layercounts_LSC(i,:) = N;
        [N,edges] = histcounts(dataglut_CSC(in_CSC,1),[1:9]);
        layercounts_CSC(i,:) = N;
        
        
    end
end

tmp = layercounts_LSC(:,1:6)./repmat(sum(layercounts_LSC(:,1:6),2),1,6)*100;
tmp(tmp==0) = 0.001;
    
figure;
set(gcf,'position',[100,100,570,1000],'color','w')
axes('position',[0.3,0.05,0.4,0.85]);
x = repmat([1:6],length(layercounts_LSC),1);
y = repmat([1:length(layercounts_LSC)]',1,6);
scatter(x(:), y(:), 6*tmp(:),'facecolor','b');
axis tight
set(gca,'xtick',[1:6],'xTickLabel',{'L1','L2','L3','L4','L5','L6'},....
    'ydir','reverse','ytick',[1:length(layercounts_LSC)],'YTickLabel',combi_needed(:,3),'xlim',[0.5,6.5],'ylim',[0.5,length(layercounts_LSC)+0.5]);
eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/blobs_clusters_vs_layers_LSC_',date,'.pdf']);

 
figure;
set(gcf,'position',[100,100,1000,600],'color','w')
axes('position',[0.05,0.3,0.9,0.3]);
x = repmat([1:6],length(layercounts_LSC),1);
y = repmat([1:length(layercounts_LSC)]',1,6);
scatter(y(:), x(:), 6*tmp(:),'facecolor','b');
axis tight
set(gca,'ytick',[1:6],'yTickLabel',{'L1','L2','L3','L4','L5','L6','CC','Out'},....
    'ydir','reverse','xtick',[1:length(layercounts_LSC)],'xTickLabel',....
    combi_needed(:,3),'ylim',[0.5,6.5],'xlim',[0.5,length(layercounts_LSC)+0.5],'XTickLabelRotation',45);
eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/blobs_clusters_vs_layers_horizon_LSC_',date,'.pdf']);
% % % % % % % % % % % % % % % % % % 

tmp = layercounts_CSC(:,1:6)./repmat(sum(layercounts_CSC(:,1:6),2),1,6)*100;
tmp(tmp==0) = 0.001;
figure;
set(gcf,'position',[100,100,570,1000],'color','w')
axes('position',[0.3,0.05,0.4,0.85]);
x = repmat([1:6],length(layercounts_CSC),1);
y = repmat([1:length(layercounts_CSC)]',1,6);
scatter(x(:), y(:), 6*tmp(:),'facecolor','b');
axis tight
set(gca,'xtick',[1:6],'xTickLabel',{'L1','L2','L3','L4','L5','L6','CC','Out'},....
    'ydir','reverse','ytick',[1:length(layercounts_CSC)],'YTickLabel',combi_needed(:,3),'xlim',[0.5,6.5],'ylim',[0.5,length(layercounts_CSC)+0.5]);
eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/blobs_clusters_vs_layers_CSC_',date,'.pdf']);

 
figure;
set(gcf,'position',[100,100,1000,600],'color','w')
axes('position',[0.05,0.3,0.9,0.3]);
x = repmat([1:6],length(layercounts_CSC),1);
y = repmat([1:length(layercounts_CSC)]',1,6);
scatter(y(:), x(:), 6*tmp(:),'facecolor','b');
axis tight
set(gca,'ytick',[1:6],'yTickLabel',{'L1','L2','L3','L4','L5','L6','CC','Out'},....
    'ydir','reverse','xtick',[1:length(layercounts_CSC)],'xTickLabel',....
    combi_needed(:,3),'ylim',[0.5,6.5],'xlim',[0.5,length(layercounts_CSC)+0.5],'XTickLabelRotation',45);
eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/blobs_clusters_vs_layers_horizon_CSC_',date,'.pdf']);
% % % % % % % % % % % % % % % % % % 




toc