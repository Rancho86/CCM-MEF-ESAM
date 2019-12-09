close all
clc
clear
[filename,pathname,q]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.png'},'select file', 'MultiSelect', 'on');  %可以选取多个图像
 n=length(filename);
for p=1:n
    sort1D=importdata(fullfile(pathname,filename{p}));
    sort2D=sort1D(:,:);
   x=sort2D(:,1);
   y=sort2D(:,2);
  R(1,:)=x;
  R(2,:)=y;
  R(3,:)=1;
  x_mean1=x(2)-x(1);
  y_mean1=y(2)-y(1);
  CHANGBIAN=sqrt(x_mean1 .^2+y_mean1 .^2);
  M=[x_mean1/CHANGBIAN y_mean1/CHANGBIAN 0;
      -y_mean1/CHANGBIAN x_mean1/CHANGBIAN 0;
      0 0 1];
  R1=M*R;
  figure;
 plot(R1(1,1),R1(2,1),'r+');
 hold on
 plot(R1(1,:),R1(2,:),'ro');
    x_mean = mean(R1(1,:));
    y_mean = mean(R1(2,:));
    x_v = R1(1,:) - x_mean;
    y_v = R1(2,:) - y_mean;
    x_junzhi=sqrt(sum(x_v.^2))/20;
    y_junzhi=sqrt(sum(y_v.^2))/20;
    if y_mean<0
        y_v=-y_v;
        end
if x_junzhi>y_junzhi
    beipaixufangxiang=x_v;
    paixuzhou=2;
else
    beipaixufangxiang=y_v;
    paixuzhou=1;
end
    [xfangxiang,indx] = sort(beipaixufangxiang);
    sort3D(:,:)=sort2D(indx,:);
    sort4D1=sort3D(1:4,:);
    [yfangxiang1,indy1] = sort(sort4D1(:,paixuzhou));
    sort5D1(:,:)=sort4D1(indy1,:);
    sort4D2=sort3D(5:8,:);
    [yfangxiang2,indy2] = sort(sort4D2(:,paixuzhou));
      sort5D2(:,:)=sort4D2(indy2,:);
    sort4D3=sort3D(9:12,:);
     [yfangxiang3,indy3] = sort(sort4D3(:,paixuzhou));
       sort5D3(:,:)=sort4D3(indy3,:);
    sort4D4=sort3D(13:16,:);
     [yfangxiang4,indy4] = sort(sort4D4(:,paixuzhou));
       sort5D4(:,:)=sort4D4(indy4,:);
    sort4D5=sort3D(17:20,:);
     [yfangxiang5,indy5] = sort(sort4D5(:,paixuzhou));
    sort5D5(:,:)=sort4D5(indy5,:);
 sort5D=[sort5D1;sort5D2;sort5D3;sort5D4;sort5D5];  
 figure;
 plot(sort5D(1,1),sort5D(1,2),'ro');
 hold on
 plot(sort5D(:,1),sort5D(:,2));
 filename1 = ['C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\HoloLens采集摄像头2019-11-28\排序后的二维坐标\',filename{p}(1:(end-4)),'.txt'];
    dlmwrite(filename1,sort5D,'delimiter', '\t','precision','%8.4f');
end