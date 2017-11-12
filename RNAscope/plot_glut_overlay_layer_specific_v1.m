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
colorvec = 'ckgmcbryyy'; 
colorvec = [242 235 3
    149 193 31
    0 141 54
    4 97 160
    54 169 225
    173 11 31
    244 100 28
    246 158 220
    246 158 220
    246 158 220]/255;
markervec = 'sopdosh...';
th = 3;
axes(ha(1))
for i=1:length(combi_needed(:,1))
    i
    isgaba = ~isempty(strfind(combi_needed{i,1},'GAB'));
    isglut = ~isempty(strfind(combi_needed{i,1},'Glu'));
    if isgaba
        
    elseif isglut
        tmp = strsplit(combi_needed{i,2},',');
        tmpfname = combi_needed{i,2};
        tmpfname = regexprep(tmpfname,'+','pos_');
        tmpfname = regexprep(tmpfname,'-','neg_');
        tmpfname = regexprep(tmpfname,',','-');
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
        rng(111);
        in_LSC = find(in_LSC);
        in_LSC = in_LSC(randperm(length(in_LSC)));
        in_LSC = in_LSC(1:min([100,length(in_LSC)]));
        if strcmpi(markervec(i),'.')
%             plot(dataglut_LSC(in_LSC,2), dataglut_LSC(in_LSC,3),markervec(i),'color',colorvec(i,:),'markersize',8); hold on;
            plot(dataglut_LSC(in_LSC,2), dataglut_LSC(in_LSC,3),[markervec(i),colorvec(i)],'color',colorvec(i,:),'markersize',27,'markerfacecolor',colorvec(i,:)); hold on;
        else
            plot(dataglut_LSC(in_LSC,2), dataglut_LSC(in_LSC,3),[markervec(i)],'color',colorvec(i,:),'markersize',9,'markerfacecolor',colorvec(i,:)); hold on;
        end
%         title([combi_needed{i,3},', ',' (',num2str(sum(in_LSC)),' cells)'],'fontsize',5);
    end
end
for ii=1:8
    for j=1:length(B_LSC{ii})
        plot(B_LSC{ii}{j}(:,2),B_LSC{ii}{j}(:,1),'color',[0.7,0.7,0.7]); hold on
    end
end
set(gca,'ydir','reverse','fontsize',12,'color',0.8*[1,1,1]);
axis tight
axis equal
axis off
hleg = legend(combi_needed(:,3));
set(hleg,'position',[0.37,0.8,0.1,0.15])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
axes(ha(2))
for i=1:length(combi_needed(:,1))
    i
    isgaba = ~isempty(strfind(combi_needed{i,1},'GAB'));
    isglut = ~isempty(strfind(combi_needed{i,1},'Glu'));
    if isgaba
        
    elseif isglut
        tmp = strsplit(combi_needed{i,2},',');
        tmpfname = combi_needed{i,2};
        tmpfname = regexprep(tmpfname,'+','pos_');
        tmpfname = regexprep(tmpfname,'-','neg_');
        tmpfname = regexprep(tmpfname,',','-');
        in_CSC = true(length(dataglut_CSC(:,1)),1);
        in_CSC = true(length(dataglut_CSC(:,1)),1);
        in_CSC(cellfun(@isempty, strfind(imagesource_glut_CSC, combi_needed{i,1}))) = false;
        in_CSC(cellfun(@isempty, strfind(imagesource_glut_CSC, combi_needed{i,1}))) = false;
        for j=1:length(tmp)
            if tmp{j}(1)=='+'
                in_CSC = in_CSC & dataglut_CSC(:,find(strcmpi(table_header_Glut,tmp{j}(2:end)) ) )>th;
                in_CSC = in_CSC & dataglut_CSC(:,find(strcmpi(table_header_Glut,tmp{j}(2:end)) ) )>th;
            elseif tmp{j}(1)=='-'
                in_CSC = in_CSC & ~(dataglut_CSC(:,find(strcmpi(table_header_Glut,tmp{j}(2:end)) ) )>th);
                in_CSC = in_CSC & ~(dataglut_CSC(:,find(strcmpi(table_header_Glut,tmp{j}(2:end)) ) )>th);
            end
        end
        rng(111);
        in_CSC = find(in_CSC);
        in_CSC = in_CSC(randperm(length(in_CSC)));
        in_CSC = in_CSC(1:min([100,length(in_CSC)]));
        if strcmpi(markervec(i),'.')
%             plot(dataglut_CSC(in_CSC,2), dataglut_CSC(in_CSC,3),markervec(i),'color',colorvec(i,:),'markersize',8); hold on;
            plot(dataglut_CSC(in_CSC,2), dataglut_CSC(in_CSC,3),[markervec(i)],'color',colorvec(i,:),'markersize',27,'markerfacecolor',colorvec(i,:)); hold on;
        else
            plot(dataglut_CSC(in_CSC,2), dataglut_CSC(in_CSC,3),[markervec(i)],'color',colorvec(i,:),'markersize',9,'markerfacecolor',colorvec(i,:)); hold on;
        end
%         title([combi_needed{i,3},', ',' (',num2str(sum(in_CSC)),' cells)'],'fontsize',5);
    end
end
for ii=1:8
    for j=1:length(B_CSC{ii})
        plot(B_CSC{ii}{j}(:,2),B_CSC{ii}{j}(:,1),'color',[0.7,0.7,0.7]); hold on
    end
end
set(gca,'ydir','reverse','fontsize',12,'color',0.8*[1,1,1]);
axis tight
axis equal
axis off
hleg = legend(combi_needed(:,3));
set(hleg,'position',[0.37,0.8,0.1,0.15])



eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/'.....
    ,'overlay_Glut_lamina_specific_',date,'.pdf']);

toc