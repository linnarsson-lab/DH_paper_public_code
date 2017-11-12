tic
clear all
close all

load /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/aggregate_arc_RNAscope_counts.mat
%dataarc_LSC dataarc_CSC imagesource_arc_LSC imagesource_arc_CSC
%table_header_arc


combi_needed = loadCellFile('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/celltypes_Arc_combinations_to_output_v2.txt');
rmv = strcmpi(combi_needed(:,3),'-') | ~cellfun(@isempty, strfind(combi_needed(:,3),'alt_'));
combi_needed(rmv,:) = [];
% combi_needed = [combi_needed(1:14,:);cell(1,4);combi_needed(15:end,:)];
% combi_needed = combi_needed(1:end-3,:);
combi_needed = combi_needed([4,3,2,1,5,6,10,7,12],:);

% [~,xi] = sort(combi_needed(:,3));
% combi_needed = combi_needed(xi,:);

load ref_LSC;
load layers_mask_LSC_ref

% figure('position',[268,86,1000,834],'color','w');
% [ha, pos] = tight_subplot(4, 4, [0.01,0.01], [0.01,0.01], [0.01,0.01]);
% ha = reshape(ha,7,4)';
% ha = ha(:);
th = 3;
th_arc = 3;
arc_perc_left_right = zeros(length(combi_needed(:,1)),2);
combi_arc_perc_left_right = zeros(length(combi_needed(:,1)),2);
for i=1:length(combi_needed(:,1))
    i
%     isarc = ~isempty(strfind(combi_needed{i,1},'Arc'));
    
    tmp = strsplit(combi_needed{i,2},',');
    tmpfname = combi_needed{i,2};
    tmpfname = regexprep(tmpfname,'+','pos_');
    tmpfname = regexprep(tmpfname,'-','neg_');
    tmpfname = regexprep(tmpfname,',','-');
    in_LSC = true(length(dataarc_LSC(:,1)),1);
    %         in_CSC = true(length(dataarc_CSC(:,1)),1);
    in_LSC(cellfun(@isempty, strfind(imagesource_arc_LSC, combi_needed{i,1}))) = false;
    %         in_CSC(cellfun(@isempty, strfind(imagesource_arc_CSC, combi_needed{i,1}))) = false;
    in_LSC_arc = in_LSC & dataarc_LSC(:,find(strcmpi(table_header_arc,'Arc') ) )>th_arc;
    
    for j=1:length(tmp)-1
        if tmp{j}(1)=='+'
            in_LSC = in_LSC & dataarc_LSC(:,find(strcmpi(table_header_arc,tmp{j}(2:end)) ) )>th;
            %                 in_CSC = in_CSC & dataarc_CSC(:,find(strcmpi(table_header_arc,tmp{j}(2:end)) ) )>th;
        elseif tmp{j}(1)=='-'
            in_LSC = in_LSC & ~(dataarc_LSC(:,find(strcmpi(table_header_arc,tmp{j}(2:end)) ) )>th);
            %                 in_CSC = in_CSC & ~(dataarc_CSC(:,find(strcmpi(table_header_arc,tmp{j}(2:end)) ) )>th);
        end
    end
    %         figure('visible','off','color','w','position',[1,1,2000,1200]);
    %         axes('position',[0.01,0.01,0.9,0.46]);
    inright = sum(dataarc_LSC(in_LSC,2)>set_ref_LSC(1,1));
    inleft = sum(dataarc_LSC(in_LSC,2)<set_ref_LSC(1,1));
    inright_arc = sum(dataarc_LSC(in_LSC_arc,2)>set_ref_LSC(1,1));
    inleft_arc = sum(dataarc_LSC(in_LSC_arc,2)<set_ref_LSC(1,1));
    
    inright_arc_overlap = sum(dataarc_LSC(in_LSC_arc & in_LSC,2)>set_ref_LSC(1,1));
    inleft_arc_overlap = sum(dataarc_LSC(in_LSC_arc & in_LSC,2)<set_ref_LSC(1,1));
    in_LSC_arc = sum(in_LSC_arc);
    
%     in_LSC = in_LSC & in_LSC_arc;
    
    arc_perc_left_right(i,:) = [100*inleft_arc_overlap/in_LSC_arc, 100*inright_arc_overlap/in_LSC_arc];
    combi_arc_perc_left_right(i,:) = [100*inleft_arc_overlap/inleft, 100*inright_arc_overlap/inright] ;
    
end

figure('position',[268,86,900,300],'color','w');
[ha, pos] = tight_subplot(2, 1, [0.01,0.01], [0.2,0.01], [0.1,0.01]);
axes(ha(2))
bar(combi_arc_perc_left_right);
axis tight;
set(gca,'xtick',[1:length(combi_needed(:,1))], 'xticklabel',combi_needed(:,3),'XTickLabelRotation',45);
ylabel('% Arc+ out of group ')
axes(ha(1))
bar(arc_perc_left_right);
axis tight;
set(gca,'xtick',[]);
ylabel('% group out of Arc')

legend('left','right')

eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/'.....
    ,'Arc_percentage_Arc_',num2str(th_arc),'_genes_th',num2str(th),'_LSC_ARC_GABAGlut_',date,'.pdf']);

figure('position',[268,86,178,300],'color','w');
[ha, pos] = tight_subplot(1, 2, [0.01,0.01], [0.1,0.01], [0.3,0.01]);
axes(ha(1))
hb = barh(combi_arc_perc_left_right); hold on;
set(hb(1),'barwidth',1.5,'facecolor','g','edgecolor','none')
set(hb(2),'barwidth',1.5,'facecolor',[255,124,36]/255,'edgecolor','none')
plot(5*[1,1],[0,length(combi_needed(:,1))+1],'color',0.6*[1,1,1],'linewidth',0.5);
plot(10*[1,1],[0,length(combi_needed(:,1))+1],'color',0.6*[1,1,1],'linewidth',0.5);
plot(15*[1,1],[0,length(combi_needed(:,1))+1],'color',0.6*[1,1,1],'linewidth',0.5);
plot(20*[1,1],[0,length(combi_needed(:,1))+1],'color',0.6*[1,1,1],'linewidth',0.5);
axis tight;
set(gca,'ytick',[1:length(combi_needed(:,1))], 'yticklabel',combi_needed(:,3),'fontsize',8,'ydir','reverse');
xlabel('% Arc+ out of group ')
box off
axes(ha(2))
hb = barh(arc_perc_left_right);
set(hb(1),'barwidth',1.5,'facecolor','g','edgecolor','none'); hold on;
set(hb(2),'barwidth',1.5,'facecolor',[255,124,36]/255,'edgecolor','none')
plot(5*[1,1],[0,length(combi_needed(:,1))+1],'color',0.6*[1,1,1],'linewidth',0.5);
plot(10*[1,1],[0,length(combi_needed(:,1))+1],'color',0.6*[1,1,1],'linewidth',0.5);
plot(15*[1,1],[0,length(combi_needed(:,1))+1],'color',0.6*[1,1,1],'linewidth',0.5);
plot(20*[1,1],[0,length(combi_needed(:,1))+1],'color',0.6*[1,1,1],'linewidth',0.5);
axis tight;
set(gca,'ytick',[],'fontsize',8,'ydir','reverse');
xlabel('% group out of Arc')
box off

legend('left','right')


eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/'.....
    ,'Arc_percentage_Arc_',num2str(th_arc),'_genes_th',num2str(th),'_LSC_ARC_GABAGlut_',date,'.pdf']);


toc