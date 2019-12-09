function f = FindRT(RTPara, A, k, M, mpic)                 
%% 用来优化外参矩阵中的 R 和 t

t = RTPara(1, 7 : 9)';
r1 = RTPara(1, 1 : 3);
r1=r1/norm(r1);
r2 = RTPara(1, 4 : 6); 
r2=r2/norm(r2);
r3 = cross(r1, r2);
R = [r1; r2; r3];
k1=k(1); k2=k(2); k3=k(3);
p1=k(4); p2=k(5);
% fu = A(1, 1); fv = A(2, 2); s = A(1, 2); 
u0 = A(1,3); v0 = A(2,3);
[m, n, numFiles] = size(M);
D = []; d = [];
for i = 1 : numFiles
    for j = 1 : n

        u = mpic(1, j, i) ;
        v =  mpic(2, j, i) ;
        %构建 DLT 矩阵
        d = [d; u; v];
        % 求理想坐标
        Xw = M(1, j, i); Yw = M(2, j, i); Zw = M(3, j, i);
        UV = A * [R, t] * [Xw; Yw; Zw; 1];
        U = UV(1) / UV(3); V = UV(2) / UV(3);
        D = [D; U; V];
    end
end
f = [D-d;( r1 * r1' - 1) +  ( r2 * r2' - 1) + r1 * r2'];%
