%Jin Z, Gupta V, Murray R M. State estimation over packet dropping networks using multiple description coding %☆[J]. Automatica, 2006, 42(9):1441-1452.
%By xuelang-wang

% function Md_code
close all
clear
clc
% A=[1.25 0;1 1.1];C=[1 1];Q=20*eye(2,2);R=2.5;       % 参数设置
A=-1.25;C=1;Q=1;R=2.5;
%%%%%%三链路
% A=[1.1 0.2 zeros(1,6);0.2 0.3 0.3 zeros(1,5);0 0.3 0.5 0.2 zeros(1,4);
%     0 0 0.2 1.1 0.3 zeros(1,3);zeros(1,3) 0.3 0.2 0.4 0 0;zeros(1,4) 0.4 0.4 0.3 0;
%     zeros(1,5) 0.3 0.5 0.3;zeros(1,6) 0.3 1.1];
% Q=20*eye(8);
% C1=[1 1 1 zeros(1,5);0 0 1 1 zeros(1,4)];R1=2.5*eye(2);
% C2=[0 0 1 1 0 0 0 0;0 0 0 0 1 1 0 0];R2=2.5*eye(2);
% C3=[zeros(1,4) 1 1 0 0;zeros(1,5) 1 1 1];R3=2.5*eye(2);
% C=[C1;C2;C3];k=3;R=2.5*eye(6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambdadown=1-1/(max(eig(A)))^2;   %gathered下界
lambdadown2=1-1/max(abs(eig(A)));   %MDcode下界
lambda=0:0.01:1;                 %设置步长
Trs=1e10*ones(size(lambda));   %设成无穷大
Trs2=1e10*ones(size(lambda));   %设成无穷大
k=1;
for i=lambda
    if i >= lambdadown2+0.01
        A2=A*(1-i);
        s2=dlyap(A2,Q);
        Trs2(k)=trace(s2);
        if i>lambdadown
            A1=A*sqrt(1-i);
            s=dlyap(A1,Q);
            Trs(k)=trace(s);
        end
    end
    k=k+1;
end
plot(lambda,Trs,'-b',lambda,Trs2,'--r')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%二分法求gathered上界
u=0;v=1;
tol=1;
while tol >= 0.001
    lambdaup = (u+v)/2;
    setlmis([]);
    Y=lmivar(1,[1,1]);       %与矩阵A同阶的对称矩阵
    Z=lmivar(2,[1,1]);       %矩阵Z的行数与矩阵A行数相等，列与矩阵C的行数相等
    lmiterm([-1 1 1 Y],1,1);
    lmiterm([-1 1 2 Y],sqrt(lambdaup),A);
    lmiterm([-1 1 2 Z],sqrt(lambdaup),C);
    lmiterm([-1 1 3 Y],sqrt(1-lambdaup),A);
    lmiterm([-1 2 2 Y],1,1);
    lmiterm([-1 3 3 Y],1,1);
    lmiterm([-2 1 1 Y],1,1);
    lmiterm([-3 1 1 0],1);
    lmiterm([3,1 1,Y],1,1);
    lmisys=getlmis;
    options=[0 0 0 0 1];
    [t,xopt]=feasp(lmisys,options);
    if t <=0
        v=lambdaup;
    else
        u=lambdaup;
    end
    tol=abs(u-v);
end
TrV=zeros(size(lambda));
k=1;
for i=lambda
    setlmis([]);
    V=lmivar(1,[1,1]);        %与矩阵A同阶的对称矩阵 
    lmiterm([-1 1 1 V],A,A');
    lmiterm([-1 1 1 V],-1,1);
    lmiterm([-1 1 1 0],Q);
    lmiterm([-1 1 2 V],sqrt(i)*A,C');
    lmiterm([-1 2 2 V],C,C');
    lmiterm([-1 2 2 0],R);
    lmiterm([-2 1 1 V],1,1);
    lmisys=getlmis;
    n=decnbr(lmisys);
    c=zeros(n,1);
    for j=1:n
        v=defcx(lmisys,j,V);
        c(j)=-1*trace(v);
    end
    options=[0 0 0 0 1];
    [copt,xopt]=mincx(lmisys,c,options);
    TrV(k)=-copt;
    k=k+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%二分法求MD code上界
u2=0;v2=1;
tol2=1;
while tol2 >= 0.001
    lambdaup2 = (u2+v2)/2;
    setlmis([]);
    Y=lmivar(1,[1,1]);       %与矩阵A同阶的对称矩阵
    Z=lmivar(2,[1,1]);       %矩阵Z的行数与矩阵A行数相等，列与矩阵C的行数相等
    Z1=lmivar(2,[1,1]);      %矩阵Z的行数与矩阵A行数相等，列与矩阵C的行数相等 
    lmiterm([-1 1 1 Y],1,1);
    lmiterm([-1 1 2 Y],sqrt(2*(1-lambdaup2)*lambdaup2),A);
    lmiterm([-1 1 2 Z1],sqrt(2*(1-lambdaup2)*lambdaup2),C);
    lmiterm([-1 1 3 Y],lambdaup2,A);
    lmiterm([-1 1 3 Z],lambdaup2,C);
    lmiterm([-1 1 4 Y],1-lambdaup2,A);
    lmiterm([-1 2 2 Y],1,1);
    lmiterm([-1 3 3 Y],1,1);
    lmiterm([-1 4 4 Y],1,1);
    lmiterm([-2 1 1 Y],1,1);
    lmiterm([-3 1 1 0],1);
    lmiterm([3,1 1,Y],1,1);
    lmisys=getlmis;
    options=[0 0 0 0 1];
    [t,xopt]=feasp(lmisys,options);
    if t <=0
        v2=lambdaup2;
    else
        u2=lambdaup2;
    end
    tol2=abs(u2-v2);
end
TrV2=zeros(size(lambda));
k=1;D0=8.33e-6;D1=1.56;
G0=R+D0;
G1=R+D1;
for i=lambda
    setlmis([]);
    V=lmivar(1,[1,1]);        %与矩阵A同阶的对称矩阵 
    lmiterm([-1 1 1 V],A,A');
    lmiterm([-1 1 1 V],-1,1);
    lmiterm([-1 1 1 0],Q);
    lmiterm([-1 1 2 V],i*A,C');
    lmiterm([-1 1 3 V],sqrt(2*i*(1-i))*A,C');
    lmiterm([-1 2 2 V],C,C');
    lmiterm([-1 2 2 0],G0);
    lmiterm([-1 3 3 V],C,C');
    lmiterm([-1 3 3 0],G1);
    lmiterm([-2 1 1 V],1,1);
    lmisys=getlmis;
    n=decnbr(lmisys);
    c=zeros(n,1);
    for j=1:n
        v=defcx(lmisys,j,V);
        c(j)=-1*trace(v);
    end
    options=[0 0 0 0 1];
    [copt,xopt]=mincx(lmisys,c,options);
    TrV2(k)=-copt;
    k=k+1;
end
hold on
plot(lambda,TrV,'-',lambda,TrV2,'g-.')
legend('Lower bound without MD','Lower bound with MD','Upper bound without MD',...
    'Upper bound with MD');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=60;                                    %高度范围
n=0:h;
lambdadown=lambdadown*ones(size(n));
lambdadown2=lambdadown2*ones(size(n));
lambdaup=lambdaup*ones(size(n));
lambdaup2=lambdaup2*ones(size(n));
hold on
plot(lambdadown,n,'.k',lambdadown2,n,'xy',lambdaup,n,'.k',lambdaup2,n,'.y')
legend('Lower bound without MD','Lower bound with MD','Upper bound without MD',...
    'Upper bound with MD','Asymptote without MD','Asymptote with MD');
%text(0.3,5,['\lambda','c=',num2str(lambdadown)])
xlabel('\lambda')
ylabel('TrS/TrV')
title('General case')
text(0.35,-10000,'$\underline{\lambda}$','interpreter','latex');
text(0.47,-10000,'$\overline{\lambda}$','interpreter','latex');
axis([0,1,0,h]);                               %%注意合理范围