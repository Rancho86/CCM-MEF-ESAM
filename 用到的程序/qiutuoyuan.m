%在不预测的情况下是完美的边缘提取程序
close all
clc
clear
[filename,pathname,q]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.png'},'select file', 'MultiSelect', 'on');  %可以选取多个图像
 if ~iscell(filename)
     filename1{1}=filename;
 else
     filename1=filename;
 end
 nn=length(filename1);
 for iii=1:nn
 img=imread(fullfile(pathname,filename1{iii}));
%img  = imread('1.bmp'); %4这张图效果最好
load('KinectcameraParams.mat');
%load('HoloLenscameraParams.mat');
u=size(img,1)/2;%图片高的一半
v=size(img,2)/2;%图片长的一半
figure,imshow(img),title('裁剪前原图');
I1 = undistortImage(img, cameraParams);
figure,imshow(I1),title('裁剪前矫正后原图');
img_gray = rgb2gray(I1);
     hold on;
     [x,y,c]=ginput(1);  %手动选取目标区域
     m(1)=x;
     n(1)=y;
     plot(x,y,'r');
     k=2;
     while(c==1)
      [x1,y1,c1]=ginput(1);
      if c1==1
      m(k)=x1;
      n(k)=y1;
      plot(x,y,'r');
      line([m(k-1) m(k)],[n(k-1) n(k)]);
      k=k+1;
      c=c1;
      else
         break
      end
     end
     line([m(k-1) m(1)],[n(k-1) n(1)]);
     BW = roipoly(img_gray,m,n);  %感兴趣区域提取的二值图像
 obj=img_gray;
     env=img_gray;
     obj(BW==0)=0;  %原图中显示目标区域
%    % 
 figure, imshow(obj), title('裁剪后');
 %thresh=graythresh(obj);
 thresh =0.8;  %这个阈值能得到指向像素原点的椭圆且不会把反光的标定板区域提取出来
bw=im2bw(obj,thresh);
 figure, imshow(bw), title('二值化后');
 L = bwlabel(bw,4);
RGB = label2rgb(L);
stats = regionprops(bw, 'basic');
%绘制感兴趣区域ROI
for i=1:size(stats)
   %  rectangle('Position',[stats(i).BoundingBox],'LineWidth',2,'LineStyle','--','EdgeColor','r'),
Xfanwei=[stats(i).BoundingBox(1)-2/3*stats(i).BoundingBox(3),stats(i).BoundingBox(1)+5/3*stats(i).BoundingBox(3),stats(i).BoundingBox(1)+5/3*stats(i).BoundingBox(3),stats(i).BoundingBox(1)-2/3*stats(i).BoundingBox(3)];
Yfanwei=[stats(i).BoundingBox(2)-2/3*stats(i).BoundingBox(4),stats(i).BoundingBox(2)-2/3*stats(i).BoundingBox(4),stats(i).BoundingBox(2)+5/3*stats(i).BoundingBox(4),stats(i).BoundingBox(2)+5/3*stats(i).BoundingBox(4)];
Xfanweimin=floor(Xfanwei(1));
Xfanweimax=floor(Xfanwei(2));
Yfanweimin=floor(Yfanwei(1));
Yfanweimax=floor(Yfanwei(3));
BW1 = roipoly(obj,Xfanwei,Yfanwei);   
 obj1=obj;
     obj1(BW1==0)=0;  %原图中显示目标区域
%     figure,imshow(obj1);
mask=[0,1,0;1,-4,1;0,1,0];%拉普拉斯滤波模板
C=imfilter(obj1,mask,'replicate');
figure, 
subplot(1,4,1),imshow(C), title('laplacian');
CC=3*C+obj1;
subplot(1,4,2), imshow(CC), title('laplacian+原图');
qq=CC(Yfanweimin:Yfanweimax,Xfanweimin:Xfanweimax);
%thresh = graythresh(qq);     %自动确定二值化阈值
thresh = 0.8;
BW2= im2bw(CC,thresh);  
subplot(1,4,3), imshow(BW2), title('laplacian+原图二值化');
imLabel = bwlabel(BW2);                %对各连通域进行标记
statsa = regionprops(imLabel,'Area');    %求各连通域的大小
area = cat(1,statsa.Area);
index = find(area == max(area));   
CCA = ismember(imLabel,index);          %获取最大连通域图像
%CCA = bwareaopen(BW2,20,8);
r=floor(1/4*(stats(i).BoundingBox(3)+stats(i).BoundingBox(4)));
%r=floor(1/5*(stats(i).BoundingBox(3)+stats(i).BoundingBox(4)));
  se=strel('disk',r);
  %fo=imopen(f,se);
 se1=strel('sphere',1);
 %  
 CCA=imfill(CCA,'holes');
 %CCA=imopen(CCA,se);
 %CCA=imclose(CCA,se1); 
%CCA=imopen(CCA,se1); 
  
 %CCA=imdilate(CCA,se1);
subplot(1,4,4), imshow(CCA), title('形态学操作');
[BBB,Q] = bwboundaries(CCA,'noholes');%寻找边缘，不包括孔
%[BBB,Q] = bwboundaries(CCA,4);%寻找边缘，不包括孔
Ba(i)=BBB;
end
figure, imshow(I1)%显示图像
tuoyuanyuanxin=[];
hold on
cuowuyuan=0;
for k = 1:length(Ba)
% x=Ba{k}
x=[];
xx=[];
   xx = Ba{k};
   for ii=1:length(xx)
   x(ii,1)=xx(ii,2);
   x(ii,2)=xx(ii,1);
   end
 p0=[0.005 0.005 0.005 0.005 0.005 0.005];
warning off
F=@(p,x)p(1)*x(:,1).^2+p(2)*x(:,1).*x(:,2)+p(3)*x(:,2).^2+p(4)*x(:,1)+p(5)*x(:,2)+p(6);
% 拟合系数，最小二乘方法
p=nlinfit(x,zeros(size(x,1),1),F,p0);
A=p(1)/p(6);
B=p(2)/p(6);
C=p(3)/p(6);
D=p(4)/p(6);
E=p(5)/p(6);
%%椭圆中心
X_center = (B*E-2*C*D)/(4*A*C - B^2);
Y_center = (B*D-2*A*E)/(4*A*C - B^2);
fprintf(' X_center=%g, Y_center=%g\n',X_center,Y_center);

%%长短轴
a= 2*sqrt((2*A*(X_center^2)+2*C*(Y_center^2)+2*B*X_center*Y_center-2)/(A+C+sqrt(((A-C)^2+B^2))));
b= 2*sqrt((2*A*(X_center^2)+2*C*(Y_center^2)+2*B*X_center*Y_center-2)/(A+C-sqrt(((A-C)^2+B^2))));
%%长轴倾角
q=0.5 * atan(B/(A-C));%这个长轴倾角可以由椭圆中心坐标和像素坐标X轴夹角求得，然后重新代入椭圆等式把B用A-C表示，重新拟合椭圆
if abs(A)<abs(C)
q=0.5 * atan(B/(A-C));
elseif abs(A)>abs(C)
   fprintf('需要加90度'); 
q=0.5 * atan(B/(A-C))+(pi/2);
end
tanzhi(k,1)=(Y_center-u)/(X_center-v);
tanzhi(k,2)=tan(q);
tanzhi(k,3)=tanzhi(k,1)-tanzhi(k,2);
i=0;
% A=0.005;
% B=0.005;
% C=0.005;
% D=0.005;
% E=0.005;
tana=2*tanzhi(k,1)/(1-tanzhi(k,1)*tanzhi(k,1));
FF=@(pp,x)pp(1)*x(:,1).^2+(pp(1)-pp(2))*tana*x(:,1).*x(:,2)+pp(2)*x(:,2).^2+pp(3)*x(:,1)+pp(4)*x(:,2)+pp(5);
% 拟合系数，最小二乘方法
 %pp0=[A B C D E];
  pp0=[0.005 0.005 0.005 0.005 0.005];
pp=nlinfit(x,zeros(size(x,1),1),FF,pp0); 
while (abs(tanzhi(k,3))>0.2)
    i=i+1;
tana=2*tanzhi(k,1)/(1-tanzhi(k,1)*tanzhi(k,1));
FF=@(pp,x)pp(1)*x(:,1).^2+(pp(1)-pp(2))*tana*x(:,1).*x(:,2)+pp(2)*x(:,2).^2+pp(3)*x(:,1)+pp(4)*x(:,2)+pp(5);
% 拟合系数，最小二乘方法
 %pp0=[A B C D E];
  pp0=[0.005 0.005 0.005 0.005 0.005];
pp=nlinfit(x,zeros(size(x,1),1),FF,pp0);
A=pp(1)/pp(5);
B=(pp(1)-pp(2))*tana/pp(5);
C=pp(2)/pp(5);
D=pp(3)/pp(5);
E=pp(4)/pp(5);
X_center = (B*E-2*C*D)/(4*A*C - B^2);
Y_center = (B*D-2*A*E)/(4*A*C - B^2);
fprintf(' 矫正后X_center=%g, Y_center=%g\n',X_center,Y_center);

%%长短轴
a= 2*sqrt((2*A*(X_center^2)+2*C*(Y_center^2)+2*B*X_center*Y_center-2)/(A+C+sqrt(((A-C)^2+B^2))));
b= 2*sqrt((2*A*(X_center^2)+2*C*(Y_center^2)+2*B*X_center*Y_center-2)/(A+C-sqrt(((A-C)^2+B^2))));
%%长轴倾角
if abs(A)<abs(C)
q=0.5 * atan(B/(A-C));
elseif abs(A)>abs(C)
   fprintf('需要加90度'); 
q=0.5 * atan(B/(A-C))+(pi/2);
end

tanzhi(k,1)=(Y_center-u)/(X_center-v);
tanzhi(k,2)=tan(q);
tanzhi(k,3)=tanzhi(k,1)-tanzhi(k,2);
if(i==51)
    cuowuyuan=cuowuyuan+1;
    break
end
end

tuoyuanyuanxin(k,1)=real(X_center);
tuoyuanyuanxin(k,2)=real(Y_center);
tuoyuanyuanxin(k,3)=real(a);
tuoyuanyuanxin(k,4)=real(b);
tuoyuanyuanxin(k,5)=real(tanzhi(k,3));
fprintf(' q=%g\n',q);
fprintf(' a=%g, b=%g\n',a,b);
plot(x(:,1),x(:,2),'r+');
hold on;
xmin=min(x(:,1));
xmax=max(x(:,1));
ymin=min(x(:,2));
ymax=max(x(:,2));
% 作图
f=ezplot(@(x,y)FF(pp,[x,y]),[xmin,xmax,ymin,ymax]);
set(f,'Color','b','LineWidth',3);
title('曲线拟合');
%legend('样本点','拟合曲线)
bb=Y_center-tanzhi(k,1)*X_center;
syms x y
[sx,sy]=solve(tanzhi(k,1)*x+bb-y==0,A*x.^2+B*x*y+C*y.^2+D*x+E*y+1==0);
if ((sx(1)-v)^2+(sy(1)-u)^2)<((sx(2)-v)^2+(sy(2)-u)^2)
tuoyuanyuanxin(k,6)=sx(1);
tuoyuanyuanxin(k,7)=sy(1);
tuoyuanyuanxin(k,8)=sx(2);
tuoyuanyuanxin(k,9)=sy(2);
elseif ((sx(1)-v)^2+(sy(1)-u)^2)>((sx(2)-v)^2+(sy(2)-u)^2)
       fprintf('需要换位置'); 
tuoyuanyuanxin(k,6)=sx(2);
tuoyuanyuanxin(k,7)=sy(2);
tuoyuanyuanxin(k,8)=sx(1);
tuoyuanyuanxin(k,9)=sy(1);
end
f=(cameraParams.FocalLength(1,1)+cameraParams.FocalLength(1,2))/2;
tanbplusc=sqrt((tuoyuanyuanxin(k,8)-v)^2+(tuoyuanyuanxin(k,9)-u)^2)/f;
tanbsubc=sqrt((tuoyuanyuanxin(k,6)-v)^2+(tuoyuanyuanxin(k,7)-u)^2)/f;
tan2b=(tanbplusc+tanbsubc)/(1-tanbplusc*tanbsubc);
beta=atan(tan2b)/2;
BILI=(tan(beta)-tanbsubc)/(tanbplusc-tanbsubc);
X_jiaozheng(k)=(tuoyuanyuanxin(k,8)-tuoyuanyuanxin(k,6))*BILI+tuoyuanyuanxin(k,6);
Y_jiaozheng(k)=(tuoyuanyuanxin(k,9)-tuoyuanyuanxin(k,7))*BILI+tuoyuanyuanxin(k,7);
qiuxinjiaozhengzuobiao(k,1)=X_jiaozheng(k);
qiuxinjiaozhengzuobiao(k,2)=Y_jiaozheng(k);
% plot(tuoyuanyuanxin(k,6),tuoyuanyuanxin(k,7),'r+')
plot([tuoyuanyuanxin(k,6) v],[tuoyuanyuanxin(k,7) u],':r','LineWidth',1);
plot([tuoyuanyuanxin(k,6) tuoyuanyuanxin(k,8)],[tuoyuanyuanxin(k,7) tuoyuanyuanxin(k,9)],'-g','LineWidth',1);
% plot(tuoyuanyuanxin(k,8),tuoyuanyuanxin(k,9),'ro')
 plot(v,u,'b+');
% plot(X_jiaozheng(k),Y_jiaozheng(k),'b.');
fprintf(' 错误圆的个数=%g\n',cuowuyuan);
end%整个循环表示的是描边
%%
%排序二维点,思路是旋转到正方向，然后根据方差判断哪个方向上有四排或者五排
 x=X_jiaozheng;
   y=Y_jiaozheng;
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
    sort3D(:,1)=qiuxinjiaozhengzuobiao(indx,1);
    sort3D(:,2)=qiuxinjiaozhengzuobiao(indx,2);
    sort4D1=sort3D(1:4,:);
    [yfangxiang1,indy1] = sort(sort4D1(:,paixuzhou));
    sort5D1(:,1)=sort4D1(indy1,1);
    sort5D1(:,2)=sort4D1(indy1,2);
    sort4D2=sort3D(5:8,:);
    [yfangxiang2,indy2] = sort(sort4D2(:,paixuzhou));
      sort5D2(:,1)=sort4D2(indy2,1);
    sort5D2(:,2)=sort4D2(indy2,2);
    sort4D3=sort3D(9:12,:);
     [yfangxiang3,indy3] = sort(sort4D3(:,paixuzhou));
       sort5D3(:,1)=sort4D3(indy3,1);
    sort5D3(:,2)=sort4D3(indy3,2);
    sort4D4=sort3D(13:16,:);
     [yfangxiang4,indy4] = sort(sort4D4(:,paixuzhou));
       sort5D4(:,1)=sort4D4(indy4,1);
    sort5D4(:,2)=sort4D4(indy4,2);
    sort4D5=sort3D(17:20,:);
     [yfangxiang5,indy5] = sort(sort4D5(:,paixuzhou));
       sort5D5(:,1)=sort4D5(indy5,1);
    sort5D5(:,2)=sort4D5(indy5,2);
 sort5D=[sort5D1;sort5D2;sort5D3;sort5D4;sort5D5];  
 %close all;
 pic=im2bw(img_gray,0);
 pic(:,:)=1;
for s=1:20
% pic(floor(sort5D(s,1)), floor(sort5D(s,2)))=0;
% pic(floor(sort5D(s,1))-1, floor(sort5D(s,2)))=0;
% pic(floor(sort5D(s,1))-2, floor(sort5D(s,2)))=0;
% pic(floor(sort5D(s,1))-3, floor(sort5D(s,2)))=0;
% pic(floor(sort5D(s,1))-4, floor(sort5D(s,2)))=0;
% pic(floor(sort5D(s,1))+1, floor(sort5D(s,2)))=0;
% pic(floor(sort5D(s,1))+2, floor(sort5D(s,2)))=0;
% pic(floor(sort5D(s,1))+3, floor(sort5D(s,2)))=0;
% pic(floor(sort5D(s,1))+4, floor(sort5D(s,2)))=0;
% pic(floor(sort5D(s,1)), floor(sort5D(s,2))-1)=0;
% pic(floor(sort5D(s,1)), floor(sort5D(s,2))-2)=0;
% pic(floor(sort5D(s,1)), floor(sort5D(s,2))-3)=0;
% pic(floor(sort5D(s,1)), floor(sort5D(s,2))-4)=0;
% pic(floor(sort5D(s,1)), floor(sort5D(s,2))+1)=0;
% pic(floor(sort5D(s,1)), floor(sort5D(s,2))+2)=0;
% pic(floor(sort5D(s,1)), floor(sort5D(s,2))+3)=0;
% pic(floor(sort5D(s,1)), floor(sort5D(s,2))+4)=0;
pic(floor(sort5D(s,2)), floor(sort5D(s,1)))=0;
pic(floor(sort5D(s,2))-1, floor(sort5D(s,1)))=0;
pic(floor(sort5D(s,2))-2, floor(sort5D(s,1)))=0;
pic(floor(sort5D(s,2))-3, floor(sort5D(s,1)))=0;
pic(floor(sort5D(s,2))-4, floor(sort5D(s,1)))=0;
pic(floor(sort5D(s,2))+1, floor(sort5D(s,1)))=0;
pic(floor(sort5D(s,2))+2, floor(sort5D(s,1)))=0;
pic(floor(sort5D(s,2))+3, floor(sort5D(s,1)))=0;
pic(floor(sort5D(s,2))+4, floor(sort5D(s,1)))=0;
pic(floor(sort5D(s,2)), floor(sort5D(s,1))-1)=0;
pic(floor(sort5D(s,2)), floor(sort5D(s,1))-2)=0;
pic(floor(sort5D(s,2)), floor(sort5D(s,1))-3)=0;
pic(floor(sort5D(s,2)), floor(sort5D(s,1))-4)=0;
pic(floor(sort5D(s,2)), floor(sort5D(s,1))+1)=0;
pic(floor(sort5D(s,2)), floor(sort5D(s,1))+2)=0;
pic(floor(sort5D(s,2)), floor(sort5D(s,1))+3)=0;
pic(floor(sort5D(s,2)), floor(sort5D(s,1))+4)=0;
end
figure
imshow(I1);
hold on
plot(sort5D(1,1),sort5D(1,2),'r+');
plot(sort5D(:,1),sort5D(:,2));
plot(sort5D(20,1),sort5D(20,2),'ro');
% adressString = ['C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\HoloLens采集摄像头2019-11-28\处理后的图片\' ,filename1{iii}(1:(end-4)),'.jpg'];
% imwrite(pic,adressString);
%  filename2 = ['C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\HoloLens采集摄像头2019-11-28\二维坐标\',filename1{iii}(1:(end-4)),'.txt'];
%     dlmwrite(filename2,sort5D,'delimiter', '\t','precision','%8.4f');
 end