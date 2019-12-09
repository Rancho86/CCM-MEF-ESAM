clc;
clear;
%  function RT = Calculate_RT(filesExtrin   A   k)
filesExtrin = {
    'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\1.txt'
   'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\2.txt'
       'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\3.txt'
%       'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\4.txt'
     'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\5.txt'
      'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\6.txt'
%    'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\7.txt'
    'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\8.txt'
    'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\9.txt'
%     'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\10.txt'
      'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\11.txt'
%      'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\12.txt'
%  'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\13.txt'
    'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\14.txt'
%      'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\15.txt'
%      'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\16.txt'
%   'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\17.txt'
%    'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\18.txt'
%      'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\19.txt'
%      'C:\Users\sushun\Desktop\论文发表参考文章\球形标定板矫正\实验\Kinect2019-11-28\联合坐标\20.txt'

    };     
                     %相机内参矩阵
load('KinectcameraParams.mat');  
B=cameraParams.IntrinsicMatrix;
A=B';
% A = [367.888839622251	0	0;
% 0 	367.689553984834	0;
% 207.255623751097	258.991369523093	1]'    ;           
 k = [0.0536964967866576;-0.0489008420760487;0.000748347826904272;0;0];
                                      %畸变系数
%  读入点的图像坐标mpic，世界坐标M
numFiles = size(filesExtrin, 1);
n = size (load(filesExtrin{1})', 2);

mpic = ones(3, n, numFiles);
M = ones(3, n, numFiles);

for i = 1 : numFiles
    % 读取u，v
    uv = load(filesExtrin{i});
    mpic(:, :, i) = [uv(:, 1 : 2)'; ones(1, n)]; %mpic是像素坐标的齐次坐标
    % 读取x，y，z
    XYZ = load(filesExtrin{i});
    M(:, :, i) = XYZ(:, 3 : 5)';       %M是世界坐标
end

% 内参和畸变参数
u0 = A(1,3); v0 = A(2,3);
k1=k(1); k2=k(2); k3=k(3);
p1=k(4); p2=k(5);

% 计算外参
D = []; d = [];
for i = 1 : numFiles
    for j = 1 : n
%         r_2 = ((mpic(1, j, i) - u0) * du )^2 +  ((mpic(2, j, i) - v0) * dv)^2;
%         r_4 = r_2^2;
%         r_6 = r_2^3;
        % 由畸变参数解出理想 u，v
        u = mpic(1, j, i) ;
        v = mpic(2, j, i) ;
        % 世界坐标
        Xw = M(1, j, i); Yw = M(2, j, i); Zw = M(3, j, i);
        % 把m34设为 1，利用非齐次方程求解   此处使用的是DLT算法
        D = [D; Xw, Yw, Zw, 1, 0, 0, 0, 0, -u * Xw, -u * Yw, -u * Zw;
                0, 0, 0, 0, Xw, Yw, Zw, 1, -v * Xw, -v * Yw, -v * Zw];
        d = [d; u; v];
    end
end
mtemp = (D' * D) \ D' * d;  %DLT算法求解初始值
mtemp = [mtemp; 1];
m34 = 1 / sqrt(sum(mtemp(9 : 11).^2));   % 平方和应该等于1，这样求得的是矩阵的系数
m=m34*mtemp;
m1 = m(1 : 3); m2 = m(5 : 7); m3 = m(9 : 11);
m14 = m(4); m24 = m(8);
% m = A*[R T] = [fuR11+u0R31 fuR12+u0R32 fuR13+u0R33 fuTx+u0Tz;
%                fvR21+u0R31 fvR22+u0R32 fvR23+u0R33 fvTy+v0Tz;
%                R31            R32         R33         Tz   ];
% 利用R的正交性来重新求取u0，v0，fu，fv
u0 = sum(m1 .* m3) / sum(m3 .* m3);
v0 = sum(m2 .* m3) / sum(m3 .* m3);
fu = norm(m1 - u0 * m3);
fv = norm(m2 - v0 * m3);
% 求旋转矩阵R，平移向量T
r1 = (m1 - u0 * m3) / fu;
r2 = (m2 - v0 * m3) / fv;
r3 = m3;
Tx = (m14 - u0 * m34) / fu;
Ty = (m24 - v0 * m34) / fv;
Tz = m34;
R = [r1'; r2'; r3'];
T = [Tx; Ty; Tz];
% 保证R为单位正交阵
[U, S, V] = svd(R);
R = U * V';
% 右手坐标系，如果R的秩为1，需要翻转
if (det(R) < 0)
    RT = -1 * [R, T];
else
    RT = [R, T];
end
 
% 优化R，T
 R = RT(:, 1 : 3);
 T = RT(:, 4);
 
 RTPara = [R(1, 1), R(1, 2), R(1, 3), R(2, 1), R(2, 2), R(2, 3), T(1), T(2), T(3)];
options = optimset('LargeScale','off','MaxFunEvals',10000,'MaxIter',100000,'Algorithm', 'levenberg-marquardt');% 使用lsqnonlin进行非线性最小二乘求解
RTtemp = lsqnonlin( @FindRT, RTPara ,[],[],options, A, k, M, mpic);

T1 = RTtemp(7 : 9)';
r1 = RTtemp(1 : 3);
r1=r1/norm(r1);
r2 = RTtemp(4 : 6);
r3 = cross(r1, r2);
r3=r3/norm(r3);
r2 = cross(r3, r1);
r2=r2/norm(r2);
R1 = [r1; r2; r3];
R_invert = R1';
T_invert = -R_invert * T1;
RT1 = [R1, T1];
inv([RT1;0 0 0 1]);
RT2 = [R_invert,T_invert]
jianyan;