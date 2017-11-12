tic
clear all
close all

load /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/aggregate_Glut_RNAscope_counts.mat
%dataglut_LSC dataglut_CSC imagesource_glut_LSC imagesource_glut_CSC
%table_header_Glut
load /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/aggregate_GABA_RNAscope_counts.mat
%datagaba_LSC datagaba_CSC imagesource_gaba_LSC imagesource_gaba_CSC
%table_header_GABA

combi_needed = loadCellFile('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/celltypes_combinations_to_output.txt');


load ref_LSC;
load ref_CSC;
load layers_mask_LSC_ref
load layers_mask_CSC_ref


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
        figure('visible','off','color','w','position',[1,1,2000,1200]);
        axes('position',[0.01,0.01,0.9,0.46]);
        plot(datagaba_LSC(in_LSC,2), datagaba_LSC(in_LSC,3),'sg'); hold on;
        for ii=1:8
            for j=1:length(B_LSC{ii})
                plot(B_LSC{ii}{j}(:,2),B_LSC{ii}{j}(:,1),'r');
            end
        end
        set(gca,'ydir','reverse');
        axis tight
        axis equal
        axis off
        title(['LSC, ',combi_needed{i,1},', ',combi_needed{i,2},' (',num2str(sum(in_LSC)),' cells)']);
        axes('position',[0.92,0.07,0.07,0.4]);
        [N,edges] = histcounts(datagaba_LSC(in_LSC,1),[1:9]);
        barh(edges(1:end-1),N);
        set(gca,'ytick',[1:8],'YTickLabel',{'L1','L2','L3','L4','L5','L6','CC','Out'},'ydir','reverse');
        xlabel('#cells')
        
        axes('position',[0.01,0.5,0.9,0.46]);
        plot(datagaba_CSC(in_CSC,2), datagaba_CSC(in_CSC,3),'sg'); hold on;
        for ii=1:8
            for j=1:length(B_CSC{ii})
                plot(B_CSC{ii}{j}(:,2),B_CSC{ii}{j}(:,1),'r');
            end
        end
        set(gca,'ydir','reverse');
        axis tight
        axis equal
        axis off
        title(['CSC, ',combi_needed{i,1},', ',combi_needed{i,2},' (',num2str(sum(in_CSC)),' cells)']);
        axes('position',[0.92,0.57,0.07,0.4]);
        [N,edges] = histcounts(datagaba_CSC(in_CSC,1),[1:9]);
        barh(edges(1:end-1),N);
        set(gca,'ytick',[1:8],'YTickLabel',{'L1','L2','L3','L4','L5','L6','CC','Out'},'ydir','reverse');
        xlabel('#cells')
        eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug21_2017/'.....
            ,combi_needed{i,1},'_',tmpfname,'_',date,'.pdf']);
        
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
        figure('visible','off','color','w','position',[1,1,2000,1200]);
        axes('position',[0.01,0.01,0.9,0.46]);
        plot(dataglut_LSC(in_LSC,2), dataglut_LSC(in_LSC,3),'sg'); hold on;
        for ii=1:8
            for j=1:length(B_LSC{ii})
                plot(B_LSC{ii}{j}(:,2),B_LSC{ii}{j}(:,1),'r');
            end
        end
        set(gca,'ydir','reverse');
        axis tight
        axis equal
        axis off
        title(['LSC, ',combi_needed{i,1},', ',combi_needed{i,2},' (',num2str(sum(in_LSC)),' cells)']);
        axes('position',[0.92,0.07,0.07,0.4]);
        [N,edges] = histcounts(dataglut_LSC(in_LSC,1),[1:9]);
        barh(edges(1:end-1),N);
        set(gca,'ytick',[1:8],'YTickLabel',{'L1','L2','L3','L4','L5','L6','CC','Out'},'ydir','reverse');
        xlabel('#cells')
        
        axes('position',[0.01,0.5,0.9,0.46]);
        plot(dataglut_CSC(in_CSC,2), dataglut_CSC(in_CSC,3),'sg'); hold on;
        for ii=1:8
            for j=1:length(B_CSC{ii})
                plot(B_CSC{ii}{j}(:,2),B_CSC{ii}{j}(:,1),'r');
            end
        end
        set(gca,'ydir','reverse');
        axis tight
        axis equal
        axis off
        title(['CSC, ',combi_needed{i,1},', ',combi_needed{i,2},' (',num2str(sum(in_CSC)),' cells)']);
        axes('position',[0.92,0.57,0.07,0.4]);
        [N,edges] = histcounts(dataglut_CSC(in_CSC,1),[1:9]);
        barh(edges(1:end-1),N);
        set(gca,'ytick',[1:8],'YTickLabel',{'L1','L2','L3','L4','L5','L6','CC','Out'},'ydir','reverse');
        xlabel('#cells')
        eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug21_2017/'.....
            ,combi_needed{i,1},'_',tmpfname,'_',date,'.pdf']);
        
    end
end

        
    










toc