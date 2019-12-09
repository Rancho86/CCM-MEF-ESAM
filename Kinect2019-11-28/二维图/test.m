clc,clear all  
circleParaXYR=[];  
 [filename,pathname,q]=uigetfile({'*.*';'*.bmp';'*.jpg';'*.png'},'select file', 'MultiSelect', 'on');  %����ѡȡ���ͼ��
num=15;
 if ~iscell(filename)
     filename1{1}=filename;
 else
     filename1=filename;
 end
 n=length(filename1);
 object=cell(n,1);  %����Ŀ�����ݿ�
 environment=cell(n,1);  %�����������ݿ�
 for i=1:n
     image=imread(fullfile(pathname,filename1{i}));
     a=rgb2gray(image);
     figure
     imshow(a);
     hold on;
     [x,y,c]=ginput(1);  %�ֶ�ѡȡĿ������
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
     BW = roipoly(a,m,n);  %����Ȥ������ȡ�Ķ�ֵͼ��
 
    obj=a;
    env=a;
    obj(BW==0)=0;  %ԭͼ����ʾĿ������
    env(BW~=0)=0;  %ԭͼ����ʾ��������
    %figure
    %imshow(obj);
    %figure
    %imshow(env);
    object{i}=obj;
    environment{i}=env;
    save(strcat('object',num2str(i),'.mat'),'obj');
    save(strcat('environment',num2str(i),'.mat'),'env');
    %close all;
 
 
[m,n,l] = size(obj);  


BW = edge(obj,'sobel');  
  
step_r = 0.1;  
step_angle = 0.1;  
minr = 4 ;  
maxr = 6;  
thresh = 0.6;  
  
[hough_space,hough_circle,para] = hough_circle(BW,step_r,step_angle,minr,maxr,thresh);  
figure(1),imshow(obj),title('ԭͼ')  
figure(2),imshow(BW),title('��Ե')  
figure(3),imshow(hough_circle),title('�����')  
  
circleParaXYR=para;  
  
%���  
fprintf(1,'\n---------------Բͳ��----------------\n');  
[r,c]=size(circleParaXYR);%r=size(circleParaXYR,1);  
fprintf(1,'  ����%d��Բ\n',r);%Բ�ĸ���  
fprintf(1,'  Բ��     �뾶\n');%Բ�ĸ���  
for n=1:r  
fprintf(1,'%d ��%d��%d��  %d\n',n,floor(circleParaXYR(n,1)),floor(circleParaXYR(n,2)),floor(circleParaXYR(n,3)));  
end  
  
%���Բ  
figure(4),imshow(a),title('����ͼ�е�Բ')  
hold on ;
plot(circleParaXYR(:,2), circleParaXYR(:,1), 'r+');  
 for k = 1 : size(circleParaXYR, 1)  
  t=0:0.01*pi:2*pi;  
  x=cos(t).*circleParaXYR(k,3)+circleParaXYR(k,2);y=sin(t).*circleParaXYR(k,3)+circleParaXYR(k,1);  
  plot(x,y,'r-');  
 end  
   
  
R_max=maxr;  
acu=zeros(R_max);  
stor =[];  
for j=1:R_max  
  for n=1:r  
   if  j == floor(circleParaXYR(n,3))  
       acu(j)= acu(j)+1;  
   end  
  end  
   stor=[stor;j,acu(j)];  
   %fprintf(1,'%d,%d\n',j,acu(j));  
end  
  
fprintf(1,'\n------------���Ӵ�С����Ŀͳ��---------\n');  
fprintf(1,'���Ӱ뾶�����Ӹ���\n');  
for j=1:R_max  
  if acu(j) > 0  
   fprintf(1,'%4d %8d\n',stor(j,1),stor(j,2));  
  end  
end  
  
fprintf(1,'----------------------------------------\n');  
figure(5),plot(stor(:,1),stor(:,2),'-k','LineWidth',2),title('������');  
xlabel('���Ӵ�С');  
ylabel('���Ӹ���');  
grid on;  
  
z=[0,10,20,30,40,50,60,70,80,90,11,35,25,42,48,40,20,75,88,94,23,10,20,30,40,78,60,76,84,95,58,10,20,30,40,50,60,70,80,90,100];%����z������  
Z=z(:);  
S=floor(abs(Z)*1);  
C=floor(abs(Z)*0.5);  
figure(6),scatter3(circleParaXYR(:,1),circleParaXYR(:,2),circleParaXYR(:,3)*7,'filled'),title('������ά���ӳ�');  
figure(7),
yuanxinhei=a;
BW=im2bw(yuanxinhei,0);
BW(:,:)=1;
for s=1:20
BW(circleParaXYR(s,1), circleParaXYR(s,2))=0;
BW(circleParaXYR(s,1)-1, circleParaXYR(s,2))=0;
BW(circleParaXYR(s,1)-2, circleParaXYR(s,2))=0;
BW(circleParaXYR(s,1)-3, circleParaXYR(s,2))=0;
BW(circleParaXYR(s,1)-4, circleParaXYR(s,2))=0;
BW(circleParaXYR(s,1)+1, circleParaXYR(s,2))=0;
BW(circleParaXYR(s,1)+2, circleParaXYR(s,2))=0;
BW(circleParaXYR(s,1)+3, circleParaXYR(s,2))=0;
BW(circleParaXYR(s,1)+4, circleParaXYR(s,2))=0;
BW(circleParaXYR(s,1), circleParaXYR(s,2)-1)=0;
BW(circleParaXYR(s,1), circleParaXYR(s,2)-2)=0;
BW(circleParaXYR(s,1), circleParaXYR(s,2)-3)=0;
BW(circleParaXYR(s,1), circleParaXYR(s,2)-4)=0;
BW(circleParaXYR(s,1), circleParaXYR(s,2)+1)=0;
BW(circleParaXYR(s,1), circleParaXYR(s,2)+2)=0;
BW(circleParaXYR(s,1), circleParaXYR(s,2)+3)=0;
BW(circleParaXYR(s,1), circleParaXYR(s,2)+4)=0;
end
imshow(BW);
%adressString = ['C:\Users\sushun\Desktop\�궨\���ݼ�ͼ��\2018-5-4\��������ڱ궨��ͼƬ\' sprintf('%0.4d', num) ',jpg'];
adressString = ['C:\Users\sushun\Desktop\�궨\���ݼ�ͼ��\2018-9-18Kinnect2�궨\������ͼƬ\' ,filename(1:(end-4)),'.jpg'];
imwrite(BW,adressString);
 filename1 = ['C:\Users\sushun\Desktop\�궨\���ݼ�ͼ��\2018-9-18Kinnect2�궨\��ά����\',filename(1:(end-4)),'.txt'];
 AX=circleParaXYR(:,1:2);
    dlmwrite(filename1,AX,'delimiter', '\t','precision','%8.4f');




 end
close all;
imshow(BW)