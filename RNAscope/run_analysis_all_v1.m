tic
clear all
close all


fname = dir('GABA/**/*.tif');
filename_cell = {fname.name}';
path_folder_cell = {fname.folder}';
path_folder_cell = cellfun(@(x) [x,'/'], path_folder_cell,'uniformoutput',false);
a = loadCellFile('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/GABA_gene_channel.txt');
a = a(2:end,:);
enhan_prct_cell = cell2mat(a(:,7:9));
p_maxtranf_cell = cell2mat(a(:,10:12));


manual_th = 3;
fig_flag = 'off';
for i=[10]%:length(filename_cell)
    fprintf(['analyzing ',path_folder_cell{i}, filename_cell{i},'\n'])
    k = find(strcmpi(a(:,2),filename_cell{i}));
    if ~isempty(k)
        chan1name = a{k,3};
        chan2name = a{k,4};
        chan3name = a{k,5};
        chan4name = a{k,6};
        RNAscope_quant_dorsalhorn_v4(path_folder_cell{i}, filename_cell{i}, ....
            chan1name,chan2name,chan3name,chan4name,manual_th,fig_flag,enhan_prct_cell(k,:),p_maxtranf_cell(k,:));
    end
    close all
end



% % % % % % % % % % % % % % % % % % % % % % % % % 


fname = dir('Glut/**/*.tif');
filename_cell = {fname.name}';
path_folder_cell = {fname.folder}';
path_folder_cell = cellfun(@(x) [x,'/'], path_folder_cell,'uniformoutput',false);
a = loadCellFile('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/Glut_gene_channel.txt');
a = a(2:end,:);
enhan_prct_cell = cell2mat(a(:,7:9));
p_maxtranf_cell = cell2mat(a(:,10:12));

manual_th = 3;
fig_flag = 'off';
for i=[86:97]%[19,20,22,24]%:length(filename_cell)
    fprintf(['analyzing ',path_folder_cell{i}, filename_cell{i},'\n'])
    k = find(strcmpi(a(:,2),filename_cell{i}));
    if ~isempty(k)
        chan1name = a{k,3};
        chan2name = a{k,4};
        chan3name = a{k,5};
        chan4name = a{k,6};
        RNAscope_quant_dorsalhorn_v4(path_folder_cell{i}, filename_cell{i}, ....
            chan1name,chan2name,chan3name,chan4name,manual_th,fig_flag,enhan_prct_cell(k,:),p_maxtranf_cell(k,:));
    end
    close all
end





toc