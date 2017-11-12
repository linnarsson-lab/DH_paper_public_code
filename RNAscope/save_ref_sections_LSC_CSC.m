tic
clear all
close all

nbit =16;

fname = '/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/GABA/Stain_GAB3_DetB_M10_LSC3/MS_Stain_GAB3_DetB_M10_LSC3.tif';
dapi = imread(fname,'index',3);
dapi = imrotate(dapi,-90);
% dapi = dapi(4630:end,1900:end);
low_in = 100/(2^nbit-1); %double(prctile(dapi(:),1))/(2^nbit-1);%
high_in = 400/(2^nbit-1); %double(prctile(dapi(:),99))/(2^nbit-1);%
dapi = imadjust(dapi,[low_in,high_in],[0,1]);
ref_LSC = dapi;
a = loadCellFile('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/GABA/Stain_GAB3_DetB_M10_LSC3/MS_Stain_GAB3_DetB_M10_LSC3_ref_points_15-Aug-2017.txt');
set_ref_LSC = cell2mat(a);
set_ref_LSC = [set_ref_LSC(1:5,1:2);set_ref_LSC(1:5,3:4)];
im_range = [min(set_ref_LSC(:,1))-100,max(set_ref_LSC(:,1))+100;...
    min(set_ref_LSC(:,2))-100,max(set_ref_LSC(:,2))+100];
ref_LSC = ref_LSC(im_range(2,1):im_range(2,2),im_range(1,1):im_range(1,2));
set_ref_LSC = [set_ref_LSC(:,1)-im_range(1,1), set_ref_LSC(:,2)-im_range(2,1)];

figure; imshow(ref_LSC); hold on;
plot(set_ref_LSC(:,1),set_ref_LSC(:,2),'or','markersize',20);

save ref_LSC ref_LSC set_ref_LSC
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

fname = '/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/Glut/Stain_Glu9_DetA_M7_CSC4/MS_Stain_Glu9_DetA_M7_CSC4.tif';
dapi = imread(fname,'index',3);
dapi = imrotate(dapi,-90);
% dapi = dapi(4630:end,1900:end);
low_in = 100/(2^nbit-1); %double(prctile(dapi(:),1))/(2^nbit-1);%
high_in = 400/(2^nbit-1); %double(prctile(dapi(:),99))/(2^nbit-1);%
dapi = imadjust(dapi,[low_in,high_in],[0,1]);
ref_CSC = dapi;
a = loadCellFile('/mnt/sanger-data2/C1_stuff/Dorsal_horn_MH/Stainings_2017/Glut/Stain_Glu9_DetA_M7_CSC4/MS_Stain_Glu9_DetA_M7_CSC4_ref_points_17-Aug-2017.txt');
set_ref_CSC = cell2mat(a);
set_ref_CSC = [set_ref_CSC(1:5,1:2);set_ref_CSC(1:5,3:4)];
im_range = [min(set_ref_CSC(:,1))-100,max(set_ref_CSC(:,1))+100;...
    min(set_ref_CSC(:,2))-100,max(set_ref_CSC(:,2))+100];
ref_CSC = ref_CSC(im_range(2,1):im_range(2,2),im_range(1,1):im_range(1,2));
set_ref_CSC = [set_ref_CSC(:,1)-im_range(1,1), set_ref_CSC(:,2)-im_range(2,1)];

figure; imshow(ref_CSC); hold on;
plot(set_ref_CSC(:,1),set_ref_CSC(:,2),'or','markersize',20);

save ref_CSC ref_CSC set_ref_CSC


toc