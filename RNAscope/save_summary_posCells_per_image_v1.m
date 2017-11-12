tic
clear all
close all

load /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/aggregate_Glut_RNAscope_counts.mat
%dataglut_LSC dataglut_CSC imagesource_glut_LSC imagesource_glut_CSC
%table_header_Glut
load /mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/aggregate_GABA_RNAscope_counts.mat
%datagaba_LSC datagaba_CSC imagesource_gaba_LSC imagesource_gaba_CSC
%table_header_GABA

th = 3;
fname_glut_LSC = unique(imagesource_glut_LSC);
fname_glut_CSC = unique(imagesource_glut_CSC);
fname_gaba_LSC = unique(imagesource_gaba_LSC);
fname_gaba_CSC = unique(imagesource_gaba_CSC);

t1 = zeros(length(fname_glut_LSC),length(table_header_Glut(5:end)));
for i=1:length(fname_glut_LSC)
    in = strcmpi(imagesource_glut_LSC,fname_glut_LSC{i});
    t1(i,:) = sum(dataglut_LSC(in,5:end)>th);
end
t2 = zeros(length(fname_glut_CSC),length(table_header_Glut(5:end)));
for i=1:length(fname_glut_CSC)
    in = strcmpi(imagesource_glut_CSC,fname_glut_CSC{i});
    t2(i,:) = sum(dataglut_CSC(in,5:end)>th);
end
tableglut = [cellfun(@(x) x(1:end-41), [fname_glut_LSC,fname_glut_CSC],'uniformoutput',false)' , m2c([t1;t2])];
tableglut = [ [{'filename'},table_header_Glut(5:end)]; tableglut ];
saveCellFile(tableglut,['glut_summary_posCells_image_',date,'.txt'])


t1 = zeros(length(fname_gaba_LSC),length(table_header_GABA(5:end)));
for i=1:length(fname_gaba_LSC)
    in = strcmpi(imagesource_gaba_LSC,fname_gaba_LSC{i});
    t1(i,:) = sum(datagaba_LSC(in,5:end)>th);
end
t2 = zeros(length(fname_gaba_CSC),length(table_header_GABA(5:end)));
for i=1:length(fname_gaba_CSC)
    in = strcmpi(imagesource_gaba_CSC,fname_gaba_CSC{i});
    t2(i,:) = sum(datagaba_CSC(in,5:end)>th);
end
tablegaba = [cellfun(@(x) x(1:end-41), [fname_gaba_LSC,fname_gaba_CSC],'uniformoutput',false)' , m2c([t1;t2])];
tablegaba = [ [{'filename'},table_header_GABA(5:end)]; tablegaba ];
saveCellFile(tablegaba,['gaba_summary_posCells_image_',date,'.txt'])




toc
