tic
close all
clear
clc
A=[1.25 0;1 1.1];C1=[1 0];C2=[1 1];C=[C1;C2];Q=20*eye(2,2);
R=2.5*eye(2);R11=2.5*eye(1);R22=2.5*eye(1);       % 参数设置
lambda=0:0.01:1;                 %设置步长
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_c1=0:0.04:1;                                             %下界
lambdadown=max(0,1-1./(((max(abs(eig(A))))^2).*(1-lambda_c1)));   
p=plot(lambda_c1,lambdadown,'*r-');
axis([0 1 0 1]);
l=1;                              %上界
lambdaup=zeros(size(lambda_c1));
for k=lambda_c1               
    n=length(lambda);     
    T=zeros(n,1);
    for i=1:n
        setlmis([]);
        Y=lmivar(1,[2,1]);       %与矩阵A同阶的对称矩阵
        Z=lmivar(2,[2,2]);       %矩阵Z的行数与矩阵A行数相等，列与矩阵C的行数相等
        Z1=lmivar(2,[2,1]);      %矩阵Z的行数与矩阵A行数相等，列与矩阵C1的行数相等 
        Z2=lmivar(2,[2,1]);      %矩阵Z的行数与矩阵A行数相等，列与矩阵C2的行数相等
        lmiterm([-1 1 1 Y],1,1);
        lmiterm([-1 1 2 Y],sqrt(k*lambda(i)),A);
        lmiterm([-1 1 2 Z],sqrt(k*lambda(i)),C);
        lmiterm([-1 1 3 Y],sqrt((1-k)*lambda(i)),A);
        lmiterm([-1 1 3 Z2],sqrt((1-k)*lambda(i)),C2);
        lmiterm([-1 1 4 Y],sqrt((1-lambda(i))*k),A);
        lmiterm([-1 1 4 Z1],sqrt((1-lambda(i))*k),C1);
        lmiterm([-1 1 5 Y],sqrt((1-lambda(i))*(1-k)),A);
        lmiterm([-1 2 2 Y],1,1);
        lmiterm([-1 3 3 Y],1,1);
        lmiterm([-1 4 4 Y],1,1);
        lmiterm([-1 5 5 Y],1,1);
        lmiterm([-2 1 1 Y],1,1);
        lmiterm([-3 1 1 0],1);
        lmiterm([3,1 1,Y],1,1);
        lmisys=getlmis;
        options=[0 0 0 0 1];
        [t,xopt]=feasp(lmisys,options);
        T(i)=t;
    end
    for j=1:n
        if T(j)<=0
            lambdaup(l)=lambda(j-1);   %%%%%可能是（j-1或j）具体看情况
            break
        end
    end
    l=l+1;
end
hold on
plot(lambda_c1,lambdaup,'-g.');
xlabel('\lambda_{1}');
ylabel('\lambda_{2}');
legend('lowerbound','upperbound');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,y]=meshgrid(0:0.04:1);
a=sqrt((1-x).*(1-y));
TrS=1e20*ones(size(a));
for i=1:length(lambda_c1)
    for j=1:length(lambdadown)
        if max(eig(A*a(i,j)))<1
            s=dlyap(A*a(i,j),Q);
            TrS(i,j)=trace(s);
        end
    end
end
figure
mesh(x,y,TrS,'edgecolor','b');
axis([0 1 0 1 0 2.5e5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=1;
for i=0:0.04:1
    g=1;
    for j=0:0.04:1
        setlmis([]);
        V=lmivar(1,[2,1]);        %与矩阵A同阶的对称矩阵 
        lmiterm([-1 1 1 V],A,A');
        lmiterm([-1 1 1 V],-1,1);
        lmiterm([-1 1 1 0],Q);
        lmiterm([-1 1 2 V],sqrt(i*j)*A,C');
        lmiterm([-1 2 2 V],C,C');
        lmiterm([-1 2 2 0],R);
        lmiterm([-1 1 3 V],sqrt(i*(1-j))*A,C1');
        lmiterm([-1 3 3 V],C1,C1');
        lmiterm([-1 3 3 0],R11);
        lmiterm([-1 1 4 V],sqrt(j*(1-i))*A,C2');
        lmiterm([-1 4 4 V],C2,C2');
        lmiterm([-1 4 4 0],R22);
        lmiterm([-2 1 1 V],1,1);
        lmisys=getlmis;
        n=decnbr(lmisys);
        c=zeros(n,1);
        for l=1:n
            v=defcx(lmisys,l,V);
            c(l)=-1*trace(v);
        end
        options=[0 0 0 0 1];
        [copt,xopt]=mincx(lmisys,c,options);
        TrV(g,f)=-copt;
        g=g+1;
    end
    f=f+1;
end
hold on
mesh(x,y,TrV,'edgecolor','r')
axis([0 1 0 1 0 2.5e5]);
xlabel('\lambda_{1}');
ylabel('\lambda_{2}');
zlabel('TrS/TrV');
title('General case');
legend('TrS','TrV');