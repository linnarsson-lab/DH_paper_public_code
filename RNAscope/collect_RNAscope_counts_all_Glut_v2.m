tic
clear all
close all


fname = dir('Glut/**/*_position_counts_per_dapi*.txt');
filename_cell = {fname.name}';
date_file = cellfun(@datenum, {fname.date}');
path_folder_cell = {fname.folder}';
path_folder_cell = cellfun(@(x) [x,'/'], path_folder_cell,'uniformoutput',false);


unifolder = unique(path_folder_cell);
filermv = [];
for i=1:length(unifolder)
    if sum(strcmpi(path_folder_cell, unifolder{i}))>1
    in = find(strcmpi(path_folder_cell, unifolder{i}));
    [~,xi] = sort(date_file(in));
    filermv = [filermv; in(xi(1:end-1))];
    end
end
filename_cell(filermv) = [];
path_folder_cell(filermv) = [];

genetmp = [];
for i=1:length(filename_cell)
    fprintf(['reading ',path_folder_cell{i}, filename_cell{i},'\n'])
    fid = fopen([path_folder_cell{i}, filename_cell{i}]);
    tline = fgetl(fid);
    tline = strsplit(tline,'\t');
    fclose(fid);
    genetmp = [genetmp, tline(4:6)];
end
geneuni = unique(genetmp);


dataglut_CSC = zeros(1000*length(filename_cell),3+length(geneuni));
imagesource_glut_CSC = cell(1000*length(filename_cell),3+length(geneuni));
k1 = 1;
for i=1:length(filename_cell)
    fprintf(['reading ',path_folder_cell{i}, filename_cell{i},'\n'])
    isCSC = ~isempty(strfind(filename_cell{i},'CSC'));
    isLSC = ~isempty(strfind(filename_cell{i},'LSC'));
    if isCSC
        fid = fopen([path_folder_cell{i}, filename_cell{i}]);
        tline = fgetl(fid);
        tline = strsplit(tline,'\t');
        a = cell2mat(textscan(fid,repmat('%f',1,6)));
        fclose(fid);
        a(sum(a(:,4:6),2)<=3,:) = [];
        k2 = k1+length(a(:,1))-1;
        [~,loc] = ismember(tline(4:6), geneuni);
        
        dataglut_CSC(k1:k2,1:3) = a(:,1:3);
        dataglut_CSC(k1:k2,3+loc) = a(:,4:6);
        imagesource_glut_CSC(k1:k2) = repmat(filename_cell(i),k2-k1+1,1);
        k1 = k2+1;
    end
end
dataglut_CSC = dataglut_CSC(1:k2,:);
dataglut_CSC = [zeros(length(dataglut_CSC),1), dataglut_CSC];
imagesource_glut_CSC = imagesource_glut_CSC(1:k2);

load layers_mask_CSC_ref % layermask_bw_LSC B_LSC Ilayers_LSC
lininds = sub2ind(size(Ilayers_CSC), round(dataglut_CSC(:,3)), round(dataglut_CSC(:,2)));

for i=1:8
    tf = ismember(lininds,find(Ilayers_CSC(:)==i));
    dataglut_CSC(tf,1) = i;
end
    

dataglut_LSC = zeros(1000*length(filename_cell),3+length(geneuni));
imagesource_glut_LSC = cell(1000*length(filename_cell),3+length(geneuni));
k1 = 1;
for i=1:length(filename_cell)
    fprintf(['reading ',path_folder_cell{i}, filename_cell{i},'\n'])
    isCSC = ~isempty(strfind(filename_cell{i},'CSC'));
    isLSC = ~isempty(strfind(filename_cell{i},'LSC'));
    if isLSC
        fid = fopen([path_folder_cell{i}, filename_cell{i}]);
        tline = fgetl(fid);
        tline = strsplit(tline,'\t');
        a = cell2mat(textscan(fid,repmat('%f',1,6)));
        fclose(fid);
        a(sum(a(:,4:6),2)<=3,:) = [];
        k2 = k1+length(a(:,1))-1;
        [~,loc] = ismember(tline(4:6), geneuni);
        
        dataglut_LSC(k1:k2,1:3) = a(:,1:3);
        dataglut_LSC(k1:k2,3+loc) = a(:,4:6);
        imagesource_glut_LSC(k1:k2) = repmat(filename_cell(i),k2-k1+1,1);
        k1 = k2+1;
    end
end
dataglut_LSC = dataglut_LSC(1:k2,:);    
dataglut_LSC = [zeros(length(dataglut_LSC),1), dataglut_LSC];
imagesource_glut_LSC = imagesource_glut_LSC(1:k2);

load layers_mask_LSC_ref % layermask_bw_LSC B_LSC Ilayers_LSC
lininds = sub2ind(size(Ilayers_LSC), round(dataglut_LSC(:,3)), round(dataglut_LSC(:,2)));

for i=1:8
    tf = ismember(lininds,find(Ilayers_LSC(:)==i));
    dataglut_LSC(tf,1) = i;
end

table_header_Glut = [{'layer','dapi_x','dapi_y','area'},geneuni];
    

save aggregate_Glut_RNAscope_counts dataglut_LSC dataglut_CSC imagesource_glut_LSC imagesource_glut_CSC table_header_Glut

toc