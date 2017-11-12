tic
clear all
close all


fname = dir('GABA/**/*_position_counts_per_dapi*.txt');
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


datagaba_CSC = zeros(1000*length(filename_cell),3+length(geneuni));
imagesource_gaba_CSC = cell(1000*length(filename_cell),3+length(geneuni));
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
        
        datagaba_CSC(k1:k2,1:3) = a(:,1:3);
        datagaba_CSC(k1:k2,3+loc) = a(:,4:6);
        imagesource_gaba_CSC(k1:k2) = repmat(filename_cell(i),k2-k1+1,1);
        k1 = k2+1;
    end
end
datagaba_CSC = datagaba_CSC(1:k2,:);
datagaba_CSC = [zeros(length(datagaba_CSC),1), datagaba_CSC];
imagesource_gaba_CSC = imagesource_gaba_CSC(1:k2);

load layers_mask_CSC_ref % layermask_bw_LSC B_LSC Ilayers_LSC
lininds = sub2ind(size(Ilayers_CSC), round(datagaba_CSC(:,3)), round(datagaba_CSC(:,2)));

for i=1:8
    tf = ismember(lininds,find(Ilayers_CSC(:)==i));
    datagaba_CSC(tf,1) = i;
end
    

datagaba_LSC = zeros(1000*length(filename_cell),3+length(geneuni));
imagesource_gaba_LSC = cell(1000*length(filename_cell),3+length(geneuni));
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
        
        datagaba_LSC(k1:k2,1:3) = a(:,1:3);
        datagaba_LSC(k1:k2,3+loc) = a(:,4:6);
        imagesource_gaba_LSC(k1:k2) = repmat(filename_cell(i),k2-k1+1,1);
        k1 = k2+1;
    end
end
datagaba_LSC = datagaba_LSC(1:k2,:);    
datagaba_LSC = [zeros(length(datagaba_LSC),1), datagaba_LSC];
imagesource_gaba_LSC = imagesource_gaba_LSC(1:k2);

load layers_mask_LSC_ref % layermask_bw_LSC B_LSC Ilayers_LSC
lininds = sub2ind(size(Ilayers_LSC), round(datagaba_LSC(:,3)), round(datagaba_LSC(:,2)));

for i=1:8
    tf = ismember(lininds,find(Ilayers_LSC(:)==i));
    datagaba_LSC(tf,1) = i;
end

table_header_GABA = [{'layer','dapi_x','dapi_y','area'},geneuni];
    

save aggregate_GABA_RNAscope_counts datagaba_LSC datagaba_CSC imagesource_gaba_LSC imagesource_gaba_CSC table_header_GABA

toc

% 
% load ref_LSC;
% in1 = find( datagaba_LSC(:,4+find(strcmpi(geneuni,'Gal')))>=0 & datagaba_LSC(:,4+find(strcmpi(geneuni,'Slc17a6')))>3....
%     & datagaba_LSC(:,4+find(strcmpi(geneuni,'Gad1')))<3 );
% in2 = find( datagaba_LSC(:,4+find(strcmpi(geneuni,'Gal')))>=0 & datagaba_LSC(:,4+find(strcmpi(geneuni,'Slc17a6')))<3....
%     & datagaba_LSC(:,4+find(strcmpi(geneuni,'Gad1')))>3 );
% in3 = find( datagaba_LSC(:,4+find(strcmpi(geneuni,'Gal')))>3 & datagaba_LSC(:,4+find(strcmpi(geneuni,'Slc17a6')))>3....
%     & datagaba_LSC(:,4+find(strcmpi(geneuni,'Gad1')))>3 );
% figure('visible','on','color','w','position',[1,1,2000,1200]);hold on;
% plot(datagaba_LSC(in1,2), datagaba_LSC(in1,3),'sg'); hold on;
% plot(datagaba_LSC(in2,2), datagaba_LSC(in2,3),'or'); hold on;
% plot(datagaba_LSC(in3,2), datagaba_LSC(in3,3),'dk'); hold on;
% plot(set_ref_LSC(:,1),set_ref_LSC(:,2),'xk','markersize',20);
% plot(set_ref_LSC(:,1),set_ref_LSC(:,2),'oc','markersize',20);
% for i=1:8
%     for j=1:length(B_LSC{i})
%         plot(B_LSC{i}{j}(:,2),B_LSC{i}{j}(:,1),'r');
%     end
% end
% set(gca,'ydir','reverse');
% axis tight
% axis equal
% axis off
% 
% 
% 
% load ref_CSC;
% in1 = find( datagaba_CSC(:,4+find(strcmpi(geneuni,'Gal')))>=0 & datagaba_CSC(:,4+find(strcmpi(geneuni,'Slc17a6')))>3....
%     & datagaba_CSC(:,4+find(strcmpi(geneuni,'Gad1')))<3 );
% in2 = find( datagaba_CSC(:,4+find(strcmpi(geneuni,'Gal')))>=0 & datagaba_CSC(:,4+find(strcmpi(geneuni,'Slc17a6')))<3....
%     & datagaba_CSC(:,4+find(strcmpi(geneuni,'Gad1')))>3 );
% in3 = find( datagaba_CSC(:,4+find(strcmpi(geneuni,'Gal')))>3 & datagaba_CSC(:,4+find(strcmpi(geneuni,'Slc17a6')))>3....
%     & datagaba_CSC(:,4+find(strcmpi(geneuni,'Gad1')))>3 );
% figure('visible','on','color','w','position',[1,1,2000,1200]);hold on;
% plot(datagaba_CSC(in1,2), datagaba_CSC(in1,3),'sg'); hold on;
% plot(datagaba_CSC(in2,2), datagaba_CSC(in2,3),'or'); hold on;
% plot(datagaba_CSC(in3,2), datagaba_CSC(in3,3),'dk'); hold on;
% plot(set_ref_CSC(:,1),set_ref_CSC(:,2),'xk','markersize',20);
% plot(set_ref_CSC(:,1),set_ref_CSC(:,2),'oc','markersize',20);
% for i=1:8
%     for j=1:length(B_CSC{i})
%         plot(B_CSC{i}{j}(:,2),B_CSC{i}{j}(:,1),'r');
%     end
% end
% set(gca,'ydir','reverse');
% axis tight
% axis equal
% axis off
% 
% 
% in1 = find( datagaba_LSC(:,4+find(strcmpi(geneuni,'Tac2')))>3 & datagaba_LSC(:,4+find(strcmpi(geneuni,'Slc17a6')))>3....
%     & datagaba_LSC(:,4+find(strcmpi(geneuni,'Gad1')))<3 );
% in2 = find( datagaba_LSC(:,4+find(strcmpi(geneuni,'Tac2')))>3 & datagaba_LSC(:,4+find(strcmpi(geneuni,'Slc17a6')))<3....
%     & datagaba_LSC(:,4+find(strcmpi(geneuni,'Gad1')))>3 );
% figure('visible','on','color','w','position',[1,1,2000,1200]);hold on;
% plot(datagaba_LSC(in1,2), datagaba_LSC(in1,3),'sg'); hold on;
% plot(datagaba_LSC(in2,2), datagaba_LSC(in2,3),'or'); hold on;
% plot(datagaba_LSC(in3,2), datagaba_LSC(in3,3),'dk'); hold on;
% plot(set_ref_LSC(:,1),set_ref_LSC(:,2),'xk','markersize',20);
% plot(set_ref_LSC(:,1),set_ref_LSC(:,2),'oc','markersize',20);
% for i=1:8
%     for j=1:length(B_LSC{i})
%         plot(B_LSC{i}{j}(:,2),B_LSC{i}{j}(:,1),'r');
%     end
% end
% set(gca,'ydir','reverse');
% axis tight
% axis equal
% axis off
% 
% 
% 
% 
