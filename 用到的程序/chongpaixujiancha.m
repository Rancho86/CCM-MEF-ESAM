close all
clear
clc
[filename,pathname,q]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.png'},'select file', 'MultiSelect', 'on');  %可以选取多个图像
 n=length(filename);
for p=1:n
sort3D=importdata(fullfile(pathname,filename{p}));
filename1 = ['C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\HoloLens采集摄像头2019-11-28\排序后的二维坐标\',filename{p}(1:(end-4)),'.txt'];
sort1D=importdata(filename1);
sort2D=[sort1D(17,:);sort1D(18,:);sort1D(19,:);sort1D(20,:);
        sort1D(13,:);sort1D(14,:);sort1D(15,:);sort1D(16,:);
        sort1D(9,:);sort1D(10,:);sort1D(11,:);sort1D(12,:);
        sort1D(5,:);sort1D(6,:);sort1D(7,:);sort1D(8,:);
        sort1D(1,:);sort1D(2,:);sort1D(3,:);sort1D(4,:)];

px=[sort2D,sort3D/1000];
    filename2 = ['C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\HoloLens采集摄像头2019-11-28\联合坐标\',filename{p}(1:(end-4)),'.txt'];
      dlmwrite(filename2,px,'delimiter', '\t','precision','%8.4f');
figure;
subplot(1,2,1);
plot(sort3D(1,3),sort3D(1,2),'ro');
hold on
plot(sort3D(:,3),sort3D(:,2));
subplot(1,2,2);
plot(sort2D(1,1),sort2D(1,2),'ro');
hold on
plot(sort2D(:,1),sort2D(:,2));
end