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

% [~,xi] = sort(combi_needed(:,3));
% combi_needed = combi_needed(xi,:);

load ref_LSC;
load layers_mask_LSC_ref

figure('position',[268,86,1350,350],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.01,0.01], [0.01,0.01], [0.01,0.01]);
% ha = reshape(ha,7,4)';
% ha = ha(:);
th_vec = [0,2,3,5,10];
for i=1:length(th_vec)
    th = th_vec(i);
    in_LSC = dataarc_LSC(:,find(strcmpi(table_header_arc,'Arc') ) )>th;
    axes(ha(i))
    for ii=1:8
        for j=1:length(B_LSC{ii})
            plot(B_LSC{ii}{j}(:,2),B_LSC{ii}{j}(:,1),'color',[0.7,0.7,0.7]);hold on;
        end
    end
    plot(dataarc_LSC(in_LSC,2), dataarc_LSC(in_LSC,3),'.r','markersize',2); hold on;
    set(gca,'ydir','reverse','fontsize',6);
    inright = dataarc_LSC(in_LSC,2)>set_ref_LSC(1,1);
    inleft = dataarc_LSC(in_LSC,2)<set_ref_LSC(1,1);
    axis tight
    axis equal
    %         set(gca,'xlim',[set_ref_LSC(1,1),set_ref_LSC(9,1)])
    axis off
    title(['Arc>',num2str(th),', ',' (',num2str(sum(in_LSC)),' cells, ',num2str(sum(inleft)),'-left, ',....
        num2str(sum(inright)),'-right[',num2str(round(100*sum(inright)/sum(in_LSC))),'%])'],'fontsize',8);
end

eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/'.....
    ,'ARC_LSC_vary_thresh_v3_',date,'.pdf']);

th_vec = [0:20];
for i=1:length(th_vec)
    th = th_vec(i);
    in_LSC = dataarc_LSC(:,find(strcmpi(table_header_arc,'Arc') ) )>th;
    
    totright(i) = sum(dataarc_LSC(in_LSC,2)>set_ref_LSC(1,1));
    totleft(i) = sum(dataarc_LSC(in_LSC,2)<set_ref_LSC(1,1));
    totalarc(i) = sum((in_LSC));
end


figure('position',[268,86,1350,350],'color','w');
[ha, pos] = tight_subplot(1, 3, [0.05,0.09], [0.09,0.03], [0.05,0.02]);
axes(ha(1))
plot(th_vec,totright,'.-');
xlabel('threshold')
ylabel('total cells or right')
axis tight
box off
axes(ha(2))
plot(th_vec,totright./totalarc,'.-');
xlabel('threshold')
ylabel('right/total')
axis tight
box off
axes(ha(3))
% plot(totright,totright./totalarc,'.-');
plot(totright/max(totright),totright./totalarc,'.-'); hold on;
d = (totright/max(totright)).^2 + (totright./totalarc).^2;
[~,imin] = min(d);
plot([0,totright(imin)/max(totright)], [0,totright(imin)./totalarc(imin)], '-r');
xlabel('total cells or right')
ylabel('right/total')
axis tight
box off
set(gca,'xlim',[0,1],'ylim',[0,1])

eval(['export_fig /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/output_celltype_distribution_aug22_2017/'.....
    ,'ARC_thresh_sensitivity_',date,'.pdf']);
