%Kothare M V, Balakrishnan V, Morai M. Robust constrained model predictive control using linear matrix %inequalities[C]// American Control Conference. IEEE, 1994:440-444 vol.1.
%example1
%By xuelang-wang

%polytope
close all
clear
clc
%参数设置
K = 0.787;%rad/(volts sec^2)
B = [0;0.1*K];%control matrix
C = [1 0];%measurement matrix
Q1 = C'*C;%量测权重
R = 0.00002;%控制权重
x0 = [0.05;0];%初始值
%%

%测试系统参数
alpha_exp = 9;
A_exp = [1 0.1;0 1-0.1*alpha_exp];
step = 40;%仿真步数
%Figure4(a)  Using norminal MPC alpha_nor = 1 sec^(-1);

%normianl system
alpha_nor = 1;%0.1 sec^(-1) <= alpha <= 10 sec(-1)
A_nor = [1 0.1;0 1-0.1*alpha_nor];
%computer control law (unconstrained infinite horizon objective function)

%use--LMItoolbox
setlmis([]);
gama = lmivar(1,[1,0]);
Q = lmivar(1,[2,1]);
Y = lmivar(2,[1,2]);

%不等式1
lmiterm([-1,1,1,0],1);
lmiterm([-1,1,2,0],x0');
lmiterm([-1,2,2,Q],1,1);

%不等式2
lmiterm([-2,1,1,Q],1,1);
lmiterm([-2,1,2,Q],1,A_nor');
lmiterm([-2,1,2,-Y],1,B');
lmiterm([-2,1,3,Q],1,sqrt(Q1));
lmiterm([-2,1,4,-Y],1,sqrt(R));
lmiterm([-2,2,2,Q],1,1);
lmiterm([-2,3,3,gama],1,1);
lmiterm([-2,4,4,gama],1,1);

%不等式3
lmiterm([-3,1,1,Q],1,1);

norminal_sys = getlmis;
c = mat2dec(norminal_sys,1,zeros(size(Q,1),size(Q,2)),zeros(size(Y,1),size(Y,2)));
[copt,xopt] = mincx(norminal_sys,c);

F_nor = dec2mat(norminal_sys,xopt,Y) / dec2mat(norminal_sys,xopt,Q)

%matlab LQR——gain
[F_dlqr,S,~] = dlqr(A_nor,B,Q1,R);
F_dlqr = -F_dlqr

%mpt toolbox
nor_model = LTISystem('A',A_nor,'B',B,'C',C);
nor_model.x.penalty = QuadFunction(Q1);
nor_model.u.penalty = QuadFunction(R);
F_mpt = nor_model.LQRGain()


x_nor = zeros(2,step+1);
x_nor(:,1) = x0;
for i = 1:40
    x_nor(:,i+1) = (A_exp+B*F_nor)*x_nor(:,i);
end
figure
plot(0:0.1:4,x_nor(1,:))
hold on;
plot(0:0.1:4,x_nor(2,:));
axis([0,4,-0.4 0.3]);
xlabel('time(sec)');
ylabel('$\theta (rad) and  \dot{\theta}(rad/sec)$','interpreter','latex');
title('Using nominal MPC with \alpha(k)=1 sec^{-1}');
h1 = legend('$\theta$','$\dot{\theta}$');
set(h1,'interpreter','latex');
%--------------------------------------------------------------------------
%Robust system
alpha_k = [0.1,10];
A_k = cell(1,2);
for i = 1:2
    A_k{1,i} = [1 0.1;0 1-0.1*alpha_k(1,i)];
end
%--------------------------------------------------------------------------
%fiugre 4(b) Using robust LMI-based MPC

x_Rot = zeros(2,step+1);
x_Rot(:,1) = x0;

    setlmis([]);
    gama = lmivar(1,[1,0]);
    Q = lmivar(1,[2,1]);
    Y = lmivar(2,[1,2]);

    %不等式1
    lmiterm([-1,1,1,0],1);
    lmiterm([-1,1,2,0],x0');
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
    Robust_sys = getlmis;
    c = mat2dec(Robust_sys,1,zeros(size(Q,1),size(Q,2)),zeros(size(Y,1),size(Y,2)));
    [copt,xopt] = mincx(Robust_sys,c);

    dec2mat(Robust_sys,xopt,gama)*inv(dec2mat(Robust_sys,xopt,Q))

    F_Rot = dec2mat(Robust_sys,xopt,Y) / dec2mat(Robust_sys,xopt,Q);

for k = 1:step
    x_Rot(:,k+1) = (A_exp + B*F_Rot)*x_Rot(:,k);
end

figure
plot(0:0.1:4,x_Rot(1,:))
hold on;
plot(0:0.1:4,x_Rot(2,:));
axis([0,4,-0.35 0.1]);
xlabel('time(sec)');
ylabel('$\theta (rad) and  \dot{\theta}(rad/sec)$','interpreter','latex');
title('Using robust LMI-based MPC');
h2 = legend('$\theta$','$\dot{\theta}$');
set(h2,'interpreter','latex');
%%
% add control constraint
umax = 2;%Euclidean norm  |u(k+i|k)| <= 2 i>=0;
step = 100;
x0 = [1;0];
x_Rot_Dynamic = zeros(2,step+1);
x_Rot_Dynamic(:,1) = x0;
u_Rot_Dynamic = zeros(1,step);

x_Rot_Static = zeros(2,step+1);
x_Rot_Static(:,1) = x0;
u_Rot_Static = zeros(1,step);

F_norm = zeros(1,step);

%computer static gain
setlmis([]);
gama = lmivar(1,[1,0]);
Q = lmivar(1,[2,1]);
Y = lmivar(2,[1,2]);

%不等式1
lmiterm([-1,1,1,0],1);
lmiterm([-1,1,2,0],x0');
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

Robust_sys_Static = getlmis;
c = mat2dec(Robust_sys_Static,1,zeros(size(Q,1),size(Q,2)),zeros(size(Y,1),size(Y,2)));
[copt,xopt] = mincx(Robust_sys_Static,c);

F_Rot_Static = dec2mat(Robust_sys_Static,xopt,Y) / dec2mat(Robust_sys_Static,xopt,Q);


for k = 1:step
    alpha_Dynamic = unifrnd(0.1,10);
    A_Dynamic = [1 0.1;0 1-0.1*alpha_Dynamic];
    
    setlmis([]);
    gama = lmivar(1,[1,0]);
    Q = lmivar(1,[2,1]);
    Y = lmivar(2,[1,2]);

    %不等式1
    lmiterm([-1,1,1,0],1);
    lmiterm([-1,1,2,0],x_Rot_Dynamic(:,k)');
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
    u_Rot_Dynamic(:,k) = F_Rot_Dynamic*x_Rot_Dynamic(:,k);
    x_Rot_Dynamic(:,k+1) = (A_Dynamic + B*F_Rot_Dynamic)*x_Rot_Dynamic(:,k);
    
    u_Rot_Static(:,k) = F_Rot_Static*x_Rot_Static(:,k);
    x_Rot_Static(:,k+1) = (A_Dynamic + B*F_Rot_Static)*x_Rot_Static(:,k);
end


figure
plot(0:0.1:step/10,x_Rot_Dynamic(1,:))
hold on
plot(0:0.1:step/10,x_Rot_Static(1,:),'--')
axis([0,10,-0.2 1]);
title('Angular position \theta (rad)');
xlabel('time(sec)');
ylabel('\theta(rad)');
legend('Dynamic','Static');

figure
plot(0:0.1:(step-1)/10,u_Rot_Dynamic(1,:));
hold on
plot(0:0.1:(step-1)/10,u_Rot_Static(1,:),'--');
axis([0,10,-2,0.5]);
title('Control signal u(volts)');
xlabel('time(sec)');
ylabel('u(volts)');
legend('Dynamic','Static');

figure
plot(0:0.1:(step-1)/10,F_norm);
hold on;
plot(0:0.1:(step-1)/10,norm(F_Rot_Static,2)*ones(100,1),'--');
% axis([0,10,0,45]);
xlabel('time(sec)');
ylabel('Norm of F');
legend('Dynamic','Static');










