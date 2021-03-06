tic
clear all
close all
load /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/aggregate_arc_RNAscope_counts.mat
%dataarc_LSC dataarc_CSC imagesource_arc_LSC imagesource_arc_CSC
%table_header_arc


combi_needed = loadCellFile('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/celltypes_combinations_to_output_v2.txt');
rmv = strcmpi(combi_needed(:,3),'-') | ~cellfun(@isempty, strfind(combi_needed(:,3),'alt_'));
combi_needed(rmv,:) = [];
combi_needed = [combi_needed(1:14,:);cell(1,4);combi_needed(15:end,:)];

% combi_needed = loadCellFile('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/celltypes_combinations_to_output_v2.txt');
% rmv = strcmpi(combi_needed(:,3),'-') | ~cellfun(@isempty, strfind(combi_needed(:,3),'alt_'));
% combi_needed(rmv,:) = [];
% combi_needed = [combi_needed(1:14,:);cell(1,4);combi_needed(15:end,:)];
% combi_needed = combi_needed([1,3,4,5,6,7,8,12,13,14],:);

combi_needed = [{'Arc8'    '-Slc17a6,+Slc32a1'        'GabaAll'};......
    {'Arc8'    '+Slc17a6,-Slc32a1'       'GlutAll'}];

% [~,xi] = sort(combi_needed(:,3));
% combi_needed = combi_needed(xi,:);

load ref_LSC;
load ref_CSC;
load layers_mask_LSC_ref
load layers_mask_CSC_ref

figure('position',[10,10,1800,950],'color','w');
[ha, pos] = tight_subplot(1, 2, [0.01,0.01], [0.01,0.01], [0.01,0.01]);
% ha = reshape(ha,7,4)';
% ha = ha(:);
% colorvec = distinguishable_colors(8);
% colorvec = [colorvec;repmat(colorvec(end,:),2,1)];
colorvec = 'rrcmgmcbryyy';
markervec = '...sopdosh...';
th = 3;

for i=1:length(combi_needed(:,1))
    i
    axes(ha(i))
%     isgaba = ~isempty(strfind(combi_needed{i,1},'GAB'));
%     isglut = ~isempty(strfind(combi_needed{i,1},'Glu'));
    
    tmp = strsplit(combi_needed{i,2},',');
    tmpfname = combi_needed{i,2};
    tmpfname = regexprep(tmpfname,'+','pos_');
    tmpfname = regexprep(tmpfname,'-','neg_');
    tmpfname = regexprep(tmpfname,',','-');
    in_LSC = true(length(dataarc_LSC(:,1)),1);
    
    in_LSC(cellfun(@isempty, strfind(imagesource_arc_LSC, 'Arc8')) ) = false;
    for j=1:length(tmp)
        if tmp{j}(1)=='+'
            in_LSC = in_LSC & dataarc_LSC(:,find(strcmpi(table_header_arc,tmp{j}(2:end)) ) )>th;
        elseif tmp{j}(1)=='-'
            in_LSC = in_LSC & ~(dataarc_LSC(:,find(strcmpi(table_header_arc,tmp{j}(2:end)) ) )>th);
        end
    end
    in_LSC = find(in_LSC);
    [bandwidth,density,X,Y] = kde2d(dataarc_LSC(in_LSC,2:3),2^6,[min(dataarc_LSC(in_LSC,2:3))-0.03*(max(dataarc_LSC(in_LSC,2:3)).....
        - min(dataarc_LSC(in_LSC,2:3)))],[max(dataarc_LSC(in_LSC,2:3))+0.03*(max(dataarc_LSC(in_LSC,2:3)) - min(dataarc_LSC(in_LSC,2:3)))]);
    contourf(X,Y,density); hold on;
    colormap('parula');

    for ii=1:8
        for j=1:length(B_LSC{ii})
            plot(B_LSC{ii}{j}(:,2),B_LSC{ii}{j}(:,1),'color',[0.7,0.7,0.7]); hold on
        end
    end
    set(gca,'ydir','reverse','fontsize',12,'color',0.8*[1,1,1]);
    axis tight
    axis equal
    axis off
end


% hleg = legend('Gad1','Slc17a6');
% set(hleg,'position',[0.37,0.8,0.1,0.15])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


eval(['export_fig -r600 /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/'.....
    ,'density_Slc32a1_Slc17a6_fromARC8_',date,'.png']);

axes(ha(1));
title('Slc32a1, LSC');
axes(ha(2));
title('Slc17a6, LSC');
eval(['export_fig -r600 /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/'.....
    ,'density_Slc32a1_Slc17a6_fromARC8_',date,'.pdf']);


% % % % % % % % % % % % % % % % % % % % % 
for i=1:length(combi_needed(:,1))
    i
    tmp = strsplit(combi_needed{i,2},',');
    tmpfname = combi_needed{i,2};
    tmpfname = regexprep(tmpfname,'+','pos_');
    tmpfname = regexprep(tmpfname,'-','neg_');
    tmpfname = regexprep(tmpfname,',','-');
    in_LSC = true(length(dataarc_LSC(:,1)),1);
    
    in_LSC(cellfun(@isempty, strfind(imagesource_arc_LSC, 'Arc8')) ) = false;
    for j=1:length(tmp)
        if tmp{j}(1)=='+'
            in_LSC = in_LSC & dataarc_LSC(:,find(strcmpi(table_header_arc,tmp{j}(2:end)) ) )>th;
        elseif tmp{j}(1)=='-'
            in_LSC = in_LSC & ~(dataarc_LSC(:,find(strcmpi(table_header_arc,tmp{j}(2:end)) ) )>th);
        end
    end
    totsum(i) = sum(in_LSC);
end

['gaba/(gaba+glut) % = ',num2str(100*totsum(1)/sum(totsum))]
['glut/(gaba+glut) % = ',num2str(100*totsum(2)/sum(totsum))]

toc