filesExtrin = {  'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\11.txt' };     
                     %相机内参矩阵
% A = [369.975614560740,0,262.099033735702;
% 0,369.364691829226,211.407670758059;
% 0,0 ,1]   ;  
load('KinectcameraParams.mat');  
B=cameraParams.IntrinsicMatrix;
A=B';
% A = [367.888839622251	0	0;
% 0 	367.689553984834	0;
% 207.255623751097	258.991369523093	1]'    ;
k = [0.0536964967866576;-0.0489008420760487;0.000748347826904272;0;0];
                                      %畸变系数
% 读入点的图像坐标mpic，世界坐标M
numFiles = size(filesExtrin, 1);
n = size (load(filesExtrin{1})', 2);

mpic = ones(3, n, numFiles);
M = ones(3, n, numFiles);

for i = 1 : numFiles
    % 读取u，v
    uv = load(filesExtrin{i});
    mpic(:, :, i) = [uv(:, 1 : 2)'; ones(1, n)];
    % 读取x，y，z
    XYZ = load(filesExtrin{i});
    M(:, :, i) = XYZ(:, 3 : 5)';
end
u1=[];
v1=[];
U1=[];
V1=[];
X1=[];
Y1=[];
Z1=[];
Xw1=[];
Yw1=[];
Zw1=[];
for i = 1 : numFiles
    for j = 1 : n
u = mpic(1, j, i) ;
        v =  mpic(2, j, i) ;
        %构建 DLT 矩阵
        %d = [d; u; v];
        % 求理想坐标
        Xw = M(1, j, i); Yw = M(2, j, i); Zw = M(3, j, i);
        UV = A * RT1 * [Xw; Yw; Zw; 1];
        U = UV(1) / UV(3); V = UV(2) / UV(3);
        %D = [D; U; V];
        u1=[u1;u];
        v1=[v1;v];
        U1=[U1;U];
        V1=[V1;V];
        
%         %对比实际坐标
%         XYZ=[u;v;1]\(A*RT1);
%         X=XYZ(1)/XYZ(4);
%         Y=XYZ(2)/XYZ(4);
%         Z=XYZ(3)/XYZ(4);
%         Xw1=[Xw1;Xw];
%         Yw1=[Yw1;Yw];
%         Zw1=[Zw1;Zw];
%         X1=[X1;X];
%         Y1=[Y1;Y];
%         Z1=[Z1;Z];
    end
end
plot(u1(:),v1(:),'ro',U1(:),V1(:),'bx')
rms=sqrt((u1(:)-U1(:)).^2+(v1(:)-V1(:)).^2);
RMS=sum(rms(:))/20 %该图上标记点标定后平均像素误差
% plot3(X1(:),Y1(:),Z1(:),'ro',Xw1(:),Yw1(:),Zw1(:),'bx')