function [K,d,P,v] = mpc_Gain_Penalty(sys,z)
%output： 
%-------u = K(x+d)----------
%K:mpc 控制增益
%d:mpc 控制常数
%------Jmin = 1/2(x-v)P(x-v)+常数
%input:
%sys:系统
%z:追踪定点

F = eye(size(sys.A));
x0 = zeros(size(sys.A,2),1);
N = 100;
[K,d,P,v]=mpc_N(sys.A,sys.B,sys.f,sys.x.penalty.H,sys.u.penalty.H,F,N,x0,z);