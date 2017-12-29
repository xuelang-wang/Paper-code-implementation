%Kothare M V, Balakrishnan V, Morai M. Robust constrained model predictive control using linear matrix %inequalities[C]// American Control Conference. IEEE, 1994:440-444 vol.1.
%example2(reference paper as following)
%B. Wie and D. S. Bernstein. Benchmark problems for robust control design. Journal of
%Guidance, Control, and Dynamics, 15(5):1057-1059, 1992
%By xuelang-wang

%polytope
close all
clear
clc
%参数设置
m1 = 1;
m2 = 1;
K_min = 0.5;%Kmin <= K <= Kmax
K_max = 10;

K_k = [0.5,10];
A_k = cell(1,2);
for i = 1:2
    A_k{1,i} = [1 0 0.1 0;0 1 0 0.1;-0.1*K_k(1,i)/m1 0.1*K_k(1,i)/m1 1 0;0.1*K_k(1,i)/m2 -0.1*K_k(1,i)/m2 0 1];
end
B = [0;0;0.1/m1;0];
C = [0 1 0 0];

Q1 = C'*C;%量测权重
R = 0.00002;%控制权重

%



% add control constraint
umax = 1;%Euclidean norm  |u(k+i|k)| <= 2 i>=0;
xs = [1;1;0;0];
us = 0;

step = 400;
K = K_min:0.5:K_max;
x0 = [0;0;0;0];
x_Rot_Dynamic = zeros(4,step+1,size(K,2));
u_Rot_Dynamic = zeros(step,size(K,2));
% y_Dynamic = zeros(step+1,size(K,2));
y_Dynamic = zeros(step+1,size(K,2));

F_norm = zeros(1,step);

% for K_k = 1:1
for K_k = 1:size(K,2);
%测试系统参数
K_exp = K(K_k);
A_exp = [1 0 0.1 0;0 1 0 0.1;-0.1*K_exp/m1 0.1*K_exp/m1 1 0;0.1*K_exp/m2 -0.1*K_exp/m2 0 1];
x_Rot_Dynamic(:,1,K_k) = x0;
y_Dynamic(1,K_k) = C*x_Rot_Dynamic(:,1,K_k);
for k = 1:step

    setlmis([]);
    gama = lmivar(1,[1,0]);
    Q = lmivar(1,[4,1]);
    Y = lmivar(2,[1,4]);

    %不等式1
    lmiterm([-1,1,1,0],1);
    lmiterm([-1,1,2,0],(x_Rot_Dynamic(:,k,K_k) - xs)');
    lmiterm([-1,2,2,Q],1,1);

     %不等式2
    lmiterm([-2,1,1,Q],1,1);
    %不等式3-4
    for j = 1:2
        lmiterm([-2-j,1,1,Q],1,1);
        lmiterm([-2-j,1,2,Q],1,A_k{1,j}');
        lmiterm([-2-j,1,2,-Y],1,B');
        lmiterm([-2-j,1,3,Q],1,sqrt(Q1));
        lmiterm([-2-j,1,4,-Y],1,sqrt(R));
        lmiterm([-2-j,2,2,Q],1,1);
        lmiterm([-2-j,3,3,gama],1,1);
        lmiterm([-2-j,4,4,gama],1,1);
    end
    %add control constraint
    lmiterm([-5,1,1,Q],1,1);
    lmiterm([-5,1,2,-Y],1,1);
    lmiterm([-5,2,2,0],umax*umax);
    
    Robust_sys_Dynamic = getlmis;
    c = mat2dec(Robust_sys_Dynamic,1,zeros(size(Q,1),size(Q,2)),zeros(size(Y,1),size(Y,2)));
    [copt,xopt] = mincx(Robust_sys_Dynamic,c);

    F_Rot_Dynamic = dec2mat(Robust_sys_Dynamic,xopt,Y) / dec2mat(Robust_sys_Dynamic,xopt,Q);
    F_norm(1,k) = norm(F_Rot_Dynamic,2);
    u_Rot_Dynamic(k,K_k) = F_Rot_Dynamic*(x_Rot_Dynamic(:,k,K_k)-xs);
    x_Rot_Dynamic(:,k+1,K_k) = (A_exp + B*F_Rot_Dynamic)*(x_Rot_Dynamic(:,k,K_k)-xs)+xs;
    y_Dynamic(k+1,K_k) = C*x_Rot_Dynamic(:,k+1,K_k);
end
end
figure
[n,m] = size(y_Dynamic);
X = repmat((0:0.1:40)',1,m);
Y = repmat(0.5:0.5:10,n,1);
mesh(X,Y,y_Dynamic);
figure
[nu,mu] = size(u_Rot_Dynamic);
X_u = repmat((0.1:0.1:40)',1,mu);
Y_u = repmat(0.5:0.5:10,nu,1);
mesh(X_u,Y_u,u_Rot_Dynamic);

% plot(0:0.1:step/10,x_Rot_Dynamic(1,:))
% 
% axis([0,10,-0.2 1]);
% title('Angular position \theta (rad)');
% xlabel('time(sec)');
% ylabel('\theta(rad)');
% legend('Dynamic','Static');
% 
% figure
% plot(0:0.1:(step-1)/10,u_Rot_Dynamic(1,:));
% 
% axis([0,10,-2,0.5]);
% title('Control signal u(volts)');
% xlabel('time(sec)');
% ylabel('u(volts)');
% legend('Dynamic','Static');











