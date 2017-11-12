tic
clear all
close all

load /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/aggregate_arc_RNAscope_counts.mat
%dataarc_LSC dataarc_CSC imagesource_arc_LSC imagesource_arc_CSC
%table_header_arc

th = 3;
fname_arc_LSC = unique(imagesource_arc_LSC);


t1 = zeros(length(fname_arc_LSC),length(table_header_arc(5:end)));
for i=1:length(fname_arc_LSC)
    in = strcmpi(imagesource_arc_LSC,fname_arc_LSC{i});
    t1(i,:) = sum(dataarc_LSC(in,5:end)>th);
end

tablearc = [cellfun(@(x) x(1:end-41), [fname_arc_LSC],'uniformoutput',false)' , m2c([t1])];
tablearc = [ [{'filename'},table_header_arc(5:end)]; tablearc ];
saveCellFile(tablearc,['arc_summary_posCells_image_',date,'.txt'])




toc
