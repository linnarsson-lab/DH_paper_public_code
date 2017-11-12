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
combi_needed = [combi_needed(1:14,:);cell(1,4);combi_needed(15:end,:)];
combi_needed = combi_needed([1,3,4,5,6,7,8,12,13,14],:);

combi_needed = [{'GAB1'    '-Slc17a6,+Gad1'        'GABA1-Glut-Gal'};......
    {'GAB3'    '+Slc17a6,-Gad1'       'GABA4-Tac2'}];

% [~,xi] = sort(combi_needed(:,3));
% combi_needed = combi_needed(xi,:);

load ref_LSC;
load ref_CSC;
load layers_mask_LSC_ref
load layers_mask_CSC_ref

figure('position',[10,10,1800,950],'color','w');
[ha, pos] = tight_subplot(2, 2, [0.01,0.01], [0.01,0.01], [0.01,0.01]);
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
    isgaba = ~isempty(strfind(combi_needed{i,1},'GAB'));
    isglut = ~isempty(strfind(combi_needed{i,1},'Glu'));
    
    tmp = strsplit(combi_needed{i,2},',');
    tmpfname = combi_needed{i,2};
    tmpfname = regexprep(tmpfname,'+','pos_');
    tmpfname = regexprep(tmpfname,'-','neg_');
    tmpfname = regexprep(tmpfname,',','-');
    in_LSC = true(length(datagaba_LSC(:,1)),1);
    in_CSC = true(length(datagaba_CSC(:,1)),1);
    
    in_LSC(cellfun(@isempty, strfind(imagesource_gaba_LSC, 'GAB1')) & cellfun(@isempty, strfind(imagesource_gaba_LSC, 'GAB3'))) = false;
    in_CSC(cellfun(@isempty, strfind(imagesource_gaba_CSC, 'GAB1')) & cellfun(@isempty, strfind(imagesource_gaba_CSC, 'GAB3'))) = false;
    for j=1:length(tmp)
        if tmp{j}(1)=='+'
            in_LSC = in_LSC & datagaba_LSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th;
            in_CSC = in_CSC & datagaba_CSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th;
        elseif tmp{j}(1)=='-'
            in_LSC = in_LSC & ~(datagaba_LSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th);
            in_CSC = in_CSC & ~(datagaba_CSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th);
        end
    end
    in_LSC = find(in_LSC);
    [bandwidth,density,X,Y] = kde2d(datagaba_LSC(in_LSC,2:3),2^10,[min(datagaba_LSC(in_LSC,2:3))-0.03*(max(datagaba_LSC(in_LSC,2:3)).....
        - min(datagaba_LSC(in_LSC,2:3)))],[max(datagaba_LSC(in_LSC,2:3))+0.03*(max(datagaba_LSC(in_LSC,2:3)) - min(datagaba_LSC(in_LSC,2:3)))]);
    contourf(X,Y,density); hold on;
    colormap('parula');
%     colormap(flipud(colormap))
%     if strcmpi(markervec(i),'.')
%         %             plot(datagaba_LSC(in_LSC,2), datagaba_LSC(in_LSC,3),markervec(i),'color',colorvec(i,:),'markersize',8); hold on;
%         plot(datagaba_LSC(in_LSC,2), datagaba_LSC(in_LSC,3),[markervec(i),colorvec(i)],'markersize',6,'markerfacecolor',colorvec(i)); hold on;
%     else
%         plot(datagaba_LSC(in_LSC,2), datagaba_LSC(in_LSC,3),[markervec(i),colorvec(i)],'markersize',3,'markerfacecolor',colorvec(i)); hold on;
%     end
    %         title([combi_needed{i,3},', ',' (',num2str(sum(in_LSC)),' cells)'],'fontsize',5);
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

% % % % % % % % % % % % % % % % % % % % % 
for i=1:length(combi_needed(:,1))
    i
    tmp = strsplit(combi_needed{i,2},',');
    tmpfname = combi_needed{i,2};
    tmpfname = regexprep(tmpfname,'+','pos_');
    tmpfname = regexprep(tmpfname,'-','neg_');
    tmpfname = regexprep(tmpfname,',','-');
    in_LSC = true(length(datagaba_LSC(:,1)),1);
       
    in_LSC(cellfun(@isempty, strfind(imagesource_gaba_LSC, 'GAB1')) & cellfun(@isempty, strfind(imagesource_gaba_LSC, 'GAB3'))) = false;
    for j=1:length(tmp)
        if tmp{j}(1)=='+'
            in_LSC = in_LSC & datagaba_LSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th;
        elseif tmp{j}(1)=='-'
            in_LSC = in_LSC & ~(datagaba_LSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th);
        end
    end
    totsum(i) = sum(in_LSC);
end

['gaba/(gaba+glut) % = ',num2str(100*totsum(1)/sum(totsum))]
['glut/(gaba+glut) % = ',num2str(100*totsum(2)/sum(totsum))]

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


for i=1:length(combi_needed(:,1))
    i
    axes(ha(i+2))
    isgaba = ~isempty(strfind(combi_needed{i,1},'GAB'));
    isglut = ~isempty(strfind(combi_needed{i,1},'Glu'));
    
    tmp = strsplit(combi_needed{i,2},',');
    tmpfname = combi_needed{i,2};
    tmpfname = regexprep(tmpfname,'+','pos_');
    tmpfname = regexprep(tmpfname,'-','neg_');
    tmpfname = regexprep(tmpfname,',','-');
    in_LSC = true(length(datagaba_LSC(:,1)),1);
    in_CSC = true(length(datagaba_CSC(:,1)),1);
    
    in_LSC(cellfun(@isempty, strfind(imagesource_gaba_LSC, 'GAB1')) & cellfun(@isempty, strfind(imagesource_gaba_LSC, 'GAB3'))) = false;
    in_CSC(cellfun(@isempty, strfind(imagesource_gaba_CSC, 'GAB1')) & cellfun(@isempty, strfind(imagesource_gaba_CSC, 'GAB3'))) = false;
    for j=1:length(tmp)
        if tmp{j}(1)=='+'
            in_LSC = in_LSC & datagaba_LSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th;
            in_CSC = in_CSC & datagaba_CSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th;
        elseif tmp{j}(1)=='-'
            in_LSC = in_LSC & ~(datagaba_LSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th);
            in_CSC = in_CSC & ~(datagaba_CSC(:,find(strcmpi(table_header_GABA,tmp{j}(2:end)) ) )>th);
        end
    end
    in_CSC = find(in_CSC);
    [bandwidth,density,X,Y] = kde2d(datagaba_CSC(in_CSC,2:3),2^10,[min(datagaba_CSC(in_CSC,2:3))-0.03*(max(datagaba_CSC(in_CSC,2:3)).....
        - min(datagaba_CSC(in_CSC,2:3)))],[max(datagaba_CSC(in_CSC,2:3))+0.03*(max(datagaba_CSC(in_CSC,2:3)) - min(datagaba_CSC(in_CSC,2:3)))]);
    contourf(X,Y,density); hold on;
    colormap('parula');
%     colormap(flipud(colormap))
%     if strcmpi(markervec(i),'.')
%         %             plot(datagaba_CSC(in_CSC,2), datagaba_CSC(in_CSC,3),markervec(i),'color',colorvec(i,:),'markersize',8); hold on;
%         plot(datagaba_CSC(in_CSC,2), datagaba_CSC(in_CSC,3),[markervec(i),colorvec(i)],'markersize',6,'markerfacecolor',colorvec(i)); hold on;
%     else
%         plot(datagaba_CSC(in_CSC,2), datagaba_CSC(in_CSC,3),[markervec(i),colorvec(i)],'markersize',3,'markerfacecolor',colorvec(i)); hold on;
%     end
    %         title([combi_needed{i,3},', ',' (',num2str(sum(in_CSC)),' cells)'],'fontsize',5);
    for ii=1:8
        for j=1:length(B_CSC{ii})
            plot(B_CSC{ii}{j}(:,2),B_CSC{ii}{j}(:,1),'color',[0.7,0.7,0.7]); hold on
        end
    end
    set(gca,'ydir','reverse','fontsize',12,'color',0.8*[1,1,1]);
    axis tight
    axis equal
    axis off
end

% eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/'.....
%     ,'overlay_GABA_Glut_density_',date,'.pdf']);


eval(['export_fig -r600 /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/'.....
    ,'overlay_GABA_Glut_density_',date,'.png']);




toc