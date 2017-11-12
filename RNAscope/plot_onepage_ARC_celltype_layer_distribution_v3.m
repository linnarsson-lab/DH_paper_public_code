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

% [~,xi] = sort(combi_needed(:,3));
% combi_needed = combi_needed(xi,:);

load ref_LSC;
load layers_mask_LSC_ref

figure('position',[268,86,1000,834],'color','w');
[ha, pos] = tight_subplot(4, 4, [0.01,0.01], [0.01,0.01], [0.01,0.01]);
% ha = reshape(ha,7,4)';
% ha = ha(:);
th = 3;
th_arc = 3;
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
    in_LSC = in_LSC & in_LSC_arc;
    for j=1:2%length(tmp)
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
    inright = dataarc_LSC(in_LSC,2)>set_ref_LSC(1,1);
    inleft = dataarc_LSC(in_LSC,2)<set_ref_LSC(1,1);
    inright_arc = dataarc_LSC(in_LSC_arc,2)>set_ref_LSC(1,1);
    inleft_arc = dataarc_LSC(in_LSC_arc,2)<set_ref_LSC(1,1);
    axes(ha(i))
    for ii=1:8
        for j=1:length(B_LSC{ii})
            plot(B_LSC{ii}{j}(:,2),B_LSC{ii}{j}(:,1),'color',[0.7,0.7,0.7]);hold on;
        end
    end
    plot(dataarc_LSC(in_LSC,2), dataarc_LSC(in_LSC,3),'.r','markersize',4); hold on;
    set(gca,'ydir','reverse','fontsize',6);
    axis tight
    axis equal
    %         set(gca,'xlim',[set_ref_LSC(1,1),set_ref_LSC(9,1)])
    axis off
    title([combi_needed{i,3},', ',' (',num2str(sum(in_LSC)),' cells, ',...
        num2str(ceil(100*sum(inleft)/sum(in_LSC_arc))),'/',num2str(ceil(100*sum(inright)/sum(in_LSC_arc))),'% left/right)'],'fontsize',8);
    
end

eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/'.....
    ,'Arc_th',num2str(th_arc),'_genes_th',num2str(th),'_LSC_ARC_GABAGlut_v3_',date,'.pdf']);



toc