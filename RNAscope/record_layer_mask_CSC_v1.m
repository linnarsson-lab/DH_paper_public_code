tic
close all
clear all


load ref_CSC;
figure;
imshow(ref_CSC); hold on;
plot(set_ref_CSC(:,1),set_ref_CSC(:,2),'xy','markersize',20);
plot(set_ref_CSC(:,1),set_ref_CSC(:,2),'oc','markersize',20);
button = 'NO';
while strcmpi(button,'NO')
    title('Click outline','color','r');
    [xoutline, youtline, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xoutline = [xoutline;xoutline(1)];
    youtline = [youtline;youtline(1)];
    xoutline = smoothn(xoutline,5);
    youtline = smoothn(youtline,5);
    plot(xoutline, youtline,'r','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end
% % % % % % % % % % % % % % % % % % 
button = 'NO';
while strcmpi(button,'NO')
    title('Click L1 left','color','r');
    [xL1left, yL1left, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xL1left = [xL1left;xL1left(1)];
    yL1left = [yL1left;yL1left(1)];
    xL1left = smoothn(xL1left,3);
    yL1left = smoothn(yL1left,3);
    plot(xL1left, yL1left,'m','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end
button = 'NO';
while strcmpi(button,'NO')
    title('Click L1 right','color','r');
    [xL1right, yL1right, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xL1right = [xL1right;xL1right(1)];
    yL1right = [yL1right;yL1right(1)];
    xL1right = smoothn(xL1right,3);
    yL1right = smoothn(yL1right,3);
    plot(xL1right, yL1right,'m','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end
% % % % % % % % % % % % % % % % % % 
button = 'NO';
while strcmpi(button,'NO')
    title('Click L2 left','color','r');
    [xL2left, yL2left, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xL2left = [xL2left;xL2left(1)];
    yL2left = [yL2left;yL2left(1)];
    xL2left = smoothn(xL2left,3);
    yL2left = smoothn(yL2left,3);
    plot(xL2left, yL2left,'g','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end
button = 'NO';
while strcmpi(button,'NO')
    title('Click L2 right','color','r');
    [xL2right, yL2right, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xL2right = [xL2right;xL2right(1)];
    yL2right = [yL2right;yL2right(1)];
    xL2right = smoothn(xL2right,3);
    yL2right = smoothn(yL2right,3);
    plot(xL2right, yL2right,'g','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end
% % % % % % % % % % % % % % % % % % 
button = 'NO';
while strcmpi(button,'NO')
    title('Click L3 left','color','r');
    [xL3left, yL3left, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xL3left = [xL3left;xL3left(1)];
    yL3left = [yL3left;yL3left(1)];
    xL3left = smoothn(xL3left,3);
    yL3left = smoothn(yL3left,3);
    plot(xL3left, yL3left,'y','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end
button = 'NO';
while strcmpi(button,'NO')
    title('Click L3 right','color','r');
    [xL3right, yL3right, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xL3right = [xL3right;xL3right(1)];
    yL3right = [yL3right;yL3right(1)];
    xL3right = smoothn(xL3right,3);
    yL3right = smoothn(yL3right,3);
    plot(xL3right, yL3right,'y','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end
% % % % % % % % % % % % % % % % % % 
button = 'NO';
while strcmpi(button,'NO')
    title('Click L4 left','color','r');
    [xL4left, yL4left, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xL4left = [xL4left;xL4left(1)];
    yL4left = [yL4left;yL4left(1)];
    xL4left = smoothn(xL4left,3);
    yL4left = smoothn(yL4left,3);
    plot(xL4left, yL4left,'b','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end
button = 'NO';
while strcmpi(button,'NO')
    title('Click L4 right','color','r');
    [xL4right, yL4right, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xL4right = [xL4right;xL4right(1)];
    yL4right = [yL4right;yL4right(1)];
    xL4right = smoothn(xL4right,3);
    yL4right = smoothn(yL4right,3);
    plot(xL4right, yL4right,'b','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end
% % % % % % % % % % % % % % % % % % 
button = 'NO';
while strcmpi(button,'NO')
    title('Click L5 left','color','r');
    [xL5left, yL5left, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xL5left = [xL5left;xL5left(1)];
    yL5left = [yL5left;yL5left(1)];
    xL5left = smoothn(xL5left,3);
    yL5left = smoothn(yL5left,3);
    plot(xL5left, yL5left,'c','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end
button = 'NO';
while strcmpi(button,'NO')
    title('Click L5 right','color','r');
    [xL5right, yL5right, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xL5right = [xL5right;xL5right(1)];
    yL5right = [yL5right;yL5right(1)];
    xL5right = smoothn(xL5right,3);
    yL5right = smoothn(yL5right,3);
    plot(xL5right, yL5right,'c','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end
% % % % % % % % % % % % % % % % % % 
button = 'NO';
while strcmpi(button,'NO')
    title('Click L6 left','color','r');
    [xL6left, yL6left, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xL6left = [xL6left;xL6left(1)];
    yL6left = [yL6left;yL6left(1)];
    xL6left = smoothn(xL6left,3);
    yL6left = smoothn(yL6left,3);
    plot(xL6left, yL6left,'w','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end
button = 'NO';
while strcmpi(button,'NO')
    title('Click L6 right','color','r');
    [xL6right, yL6right, button, ax] = ginputc(1000,'color','w','ConnectPoints',true,'ShowPoints',true);
    xL6right = [xL6right;xL6right(1)];
    yL6right = [yL6right;yL6right(1)];
    xL6right = smoothn(xL6right,3);
    yL6right = smoothn(yL6right,3);
    plot(xL6right, yL6right,'w','linewidth',2);    
    button = questdlg('is selection ok?','test','OK','NO','defualt');
end

save CSC_Layer_mask_points xoutline youtline xL1right yL1right...
    xL2right yL2right...
    xL3right yL3right...
    xL4right yL4right...
    xL5right yL5right...
    xL6right yL6right



toc