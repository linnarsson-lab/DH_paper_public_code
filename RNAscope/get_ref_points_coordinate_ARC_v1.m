tic
clear all
close all


fname = dir('Arc_Fos/Stain*/*.tif');
filename_cell = {fname.name}';
path_folder_cell = {fname.folder}';
path_folder_cell = cellfun(@(x) [x,'/'], path_folder_cell,'uniformoutput',false);
a = loadCellFile('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/Arc_gene_channel.txt');
a = a(2:end,:);
% a = a(5:5:end,:);
enhan_prct_cell = cell2mat(a(:,7:9));
p_maxtranf_cell = cell2mat(a(:,10:12));


manual_th = 3;
fig_flag = 'on';
for i=1:length(filename_cell)
    fprintf(['analyzing ',path_folder_cell{i}, filename_cell{i},'\n'])
    k = find(strcmpi(a(:,2),filename_cell{i}));    
    if ~isempty(k) & (isempty(dir([path_folder_cell{i},'*ref_points*txt'])))
        nbit = 16;
        dapi = imread([path_folder_cell{i}, filename_cell{i}],'index',3);
        dapi = imrotate(dapi,-90);
        % dapi = dapi(4630:end,1900:end);
        low_in = 100/(2^nbit-1); %double(prctile(dapi(:),1))/(2^nbit-1);%
        high_in = 400/(2^nbit-1); %double(prctile(dapi(:),99))/(2^nbit-1);%
        dapi = imadjust(dapi,[low_in,high_in],[0,1]);
        geneC4 = imread([path_folder_cell{i}, filename_cell{i}],'index',4);
        geneC4 = imrotate(geneC4,-90);
        low_in = 100/(2^nbit-1); %double(prctile(dapi(:),1))/(2^nbit-1);%
        high_in = 140/(2^nbit-1); %double(prctile(dapi(:),99))/(2^nbit-1);%
        geneC4 = imadjust(geneC4,[low_in,high_in],[0,1]);
        dapi = 0.3*dapi + 0.7*geneC4;
        button = 'NO';
        while strcmpi(button,'NO')
            close all
            figure('color','w');
            imshow(dapi); hold on;
            text(10,500,['(',num2str(i),') ',regexprep(filename_cell{i},'_','-')],'fontsize',20,'color','w');
            title('Click left side dorsal horn (6 points)','color','r');
            [xleft, yleft, button, ax] = ginputc(6,'color','w','ConnectPoints',true,'ShowPoints',true);
            plot(xleft,yleft,'r','linewidth',2);
            title('Click right side dorsal horn (6 points)','color','g')
            [xright, yright, button, ax] = ginputc(6,'color','w','ConnectPoints',true,'ShowPoints',true);
            plot(xright,yright,'g','linewidth',2);
            button = questdlg('is selection ok?','test','OK','NO','defualt');
            table1 = m2c([xleft, yleft,xright,yright]);
            saveCellFile(table1,[path_folder_cell{i}, filename_cell{i}(1:end-4),'_ref_points_',date,'.txt']);
        end
    end
    close all
end













