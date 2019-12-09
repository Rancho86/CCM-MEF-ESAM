clc;
clear;
%  function RT = Calculate_RT(filesExtrin   A   k)
filesExtrin = {
    'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\1.txt'
   'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\2.txt'
       'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\3.txt'
%       'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\4.txt'
     'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\5.txt'
      'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\6.txt'
%    'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\7.txt'
    'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\8.txt'
    'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\9.txt'
%     'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\10.txt'
      'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\11.txt'
%      'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\12.txt'
%  'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\13.txt'
    'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\14.txt'
%      'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\15.txt'
%      'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\16.txt'
%   'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\17.txt'
%    'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\18.txt'
%      'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\19.txt'
%      'C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\Kinect2019-11-28\��������\20.txt'

    };     
                     %����ڲξ���
load('KinectcameraParams.mat');  
B=cameraParams.IntrinsicMatrix;
A=B';
% A = [367.888839622251	0	0;
% 0 	367.689553984834	0;
% 207.255623751097	258.991369523093	1]'    ;           
 k = [0.0536964967866576;-0.0489008420760487;0.000748347826904272;0;0];
                                      %����ϵ��
%  ������ͼ������mpic����������M
numFiles = size(filesExtrin, 1);
n = size (load(filesExtrin{1})', 2);

mpic = ones(3, n, numFiles);
M = ones(3, n, numFiles);

for i = 1 : numFiles
    % ��ȡu��v
    uv = load(filesExtrin{i});
    mpic(:, :, i) = [uv(:, 1 : 2)'; ones(1, n)]; %mpic������������������
    % ��ȡx��y��z
    XYZ = load(filesExtrin{i});
    M(:, :, i) = XYZ(:, 3 : 5)';       %M����������
end

% �ڲκͻ������
u0 = A(1,3); v0 = A(2,3);
k1=k(1); k2=k(2); k3=k(3);
p1=k(4); p2=k(5);

% �������
D = []; d = [];
for i = 1 : numFiles
    for j = 1 : n
%         r_2 = ((mpic(1, j, i) - u0) * du )^2 +  ((mpic(2, j, i) - v0) * dv)^2;
%         r_4 = r_2^2;
%         r_6 = r_2^3;
        % �ɻ������������� u��v
        u = mpic(1, j, i) ;
        v = mpic(2, j, i) ;
        % ��������
        Xw = M(1, j, i); Yw = M(2, j, i); Zw = M(3, j, i);
        % ��m34��Ϊ 1�����÷���η������   �˴�ʹ�õ���DLT�㷨
        D = [D; Xw, Yw, Zw, 1, 0, 0, 0, 0, -u * Xw, -u * Yw, -u * Zw;
                0, 0, 0, 0, Xw, Yw, Zw, 1, -v * Xw, -v * Yw, -v * Zw];
        d = [d; u; v];
    end
end
mtemp = (D' * D) \ D' * d;  %DLT�㷨����ʼֵ
mtemp = [mtemp; 1];
m34 = 1 / sqrt(sum(mtemp(9 : 11).^2));   % ƽ����Ӧ�õ���1��������õ��Ǿ����ϵ��
m=m34*mtemp;
m1 = m(1 : 3); m2 = m(5 : 7); m3 = m(9 : 11);
m14 = m(4); m24 = m(8);
% m = A*[R T] = [fuR11+u0R31 fuR12+u0R32 fuR13+u0R33 fuTx+u0Tz;
%                fvR21+u0R31 fvR22+u0R32 fvR23+u0R33 fvTy+v0Tz;
%                R31            R32         R33         Tz   ];
% ����R����������������ȡu0��v0��fu��fv
u0 = sum(m1 .* m3) / sum(m3 .* m3);
v0 = sum(m2 .* m3) / sum(m3 .* m3);
fu = norm(m1 - u0 * m3);
fv = norm(m2 - v0 * m3);
% ����ת����R��ƽ������T
r1 = (m1 - u0 * m3) / fu;
r2 = (m2 - v0 * m3) / fv;
r3 = m3;
Tx = (m14 - u0 * m34) / fu;
Ty = (m24 - v0 * m34) / fv;
Tz = m34;
R = [r1'; r2'; r3'];
T = [Tx; Ty; Tz];
% ��֤RΪ��λ������
[U, S, V] = svd(R);
R = U * V';
% ��������ϵ�����R����Ϊ1����Ҫ��ת
if (det(R) < 0)
    RT = -1 * [R, T];
else
    RT = [R, T];
end
 
% �Ż�R��T
 R = RT(:, 1 : 3);
 T = RT(:, 4);
 
 RTPara = [R(1, 1), R(1, 2), R(1, 3), R(2, 1), R(2, 2), R(2, 3), T(1), T(2), T(3)];
options = optimset('LargeScale','off','MaxFunEvals',10000,'MaxIter',100000,'Algorithm', 'levenberg-marquardt');% ʹ��lsqnonlin���з�������С�������
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