clear
clc
for i=2
filename = ['C:\Users\sushun\Desktop\���ķ���ο�����\���α궨�����\ʵ��\HoloLens�ɼ�����ͷ2019-11-28\��ά����\',num2str(i),'.txt'];
    ex=importdata(filename);
     A(1,1:3)=ex(1,:);
     A(1,4:6)=ex(2,:);
     A(1,7:9)=ex(3,:);
   B(i,1)=sqrt((A(i,1)- A(i,4)).^2+(A(i,2)- A(i,5)).^2+(A(i,3)- A(i,6)).^2);  
   B(i,2)= sqrt((A(i,1)- A(i,7)).^2+(A(i,2)- A(i,8)).^2+(A(i,3)- A(i,9)).^2); 
   B(i,3)= sqrt((A(i,4)- A(i,7)).^2+(A(i,5)- A(i,8)).^2+(A(i,6)- A(i,9)).^2);  
   %C(i,1:3)���м��Ǹ��㣬C(i,4:6)�����м��Ǹ�����ĵ㣬C(i,7:9)�����м��Ǹ���Զ�ĵ�
if  (78>B(i,1)&&B(i,1)>74) 
    if(121>B(i,2)&&B(i,2)>118) 
         C(i,1:3)=A(i,1:3);
         C(i,4:6)=A(i,4:6);
         C(i,7:9)=A(i,7:9);
    elseif(195>B(i,2)&&B(i,2)>191)
         C(i,1:3)=A(i,4:6);
         C(i,4:6)=A(i,1:3);
         C(i,7:9)=A(i,7:9);
    else
         fprintf('�ڶ���erro!');
    end;
elseif (121>B(i,1)&&B(i,1)>118) 
       if(78>B(i,2)&&B(i,2)>74) 
         C(i,1:3)=A(i,1:3);
         C(i,4:6)=A(i,7:9);
         C(i,7:9)=A(i,4:6);
       elseif(195>B(i,2)&&B(i,2)>191) 
         C(i,1:3)=A(i,4:6);
         C(i,4:6)=A(i,7:9);
         C(i,7:9)=A(i,1:3);
       else
         fprintf('�ڶ���erro!');
        end;
elseif (195>B(i,1)&&B(i,1)>191) 
    if(78>B(i,2)&&B(i,2)>74) 
         C(i,1:3)=A(i,7:9);
         C(i,4:6)=A(i,1:3);
         C(i,7:9)=A(i,4:6);
    elseif(121>B(i,2)&&B(i,2)>118) 
         C(i,1:3)=A(i,7:9);
         C(i,4:6)=A(i,4:6);
         C(i,7:9)=A(i,1:3);
          else
         fprintf('�ڶ���erro!');
        end;
 else
  fprintf("��һ��erro!");
      end
AA(1,1)=(C(i,4)+C(i,7))/2-C(i,1);    
AA(1,2)=(C(i,5)+C(i,8))/2-C(i,2);
AA(1,3)=(C(i,6)+C(i,9))/2-C(i,3);
BB(1,1)=C(i,4)-C(i,1);
BB(1,2)=C(i,5)-C(i,2);
BB(1,3)=C(i,6)-C(i,3);
CC=cross(AA,BB);
CC=CC/norm(CC);
AA=AA/norm(AA);
BB=cross(CC,AA);
BB=BB/norm(BB);
Rktn=[AA;BB;CC]';
Tktn=C(i,1:3)'/1000;

      TRKinectTooltoNDI(:,:,i)=[Rktn,Tktn;0 0 0 1]
   %   inv(TRKinectTooltoNDI(:,:,i));
%       TRNDItoKinectToolto(:,:,k)=inv(TRKinectTooltoNDI(:,:,k));
end