function [Kf,g,Pf,h]=mpc_N(A,B,f,Q,R,F,N,x0,z)
%有限步离散系统线性二次型状态追踪问题
%mpc控制
%x(k+1) = Ax(k)+Bu(k)+f
%J =
%1/2((x'(0)-z)Q(x'(0)-z)+u'(0)Ru(0)+...+(x'(N-1)-z)Q(x'(N-1)-z)+u'(N-1)Ru(N-1))+1/2(x'(N)-z)F(x'(N)-z)
%-------u = Kf(x+g)--------
%Jmin = 1/2(x'(0)-h)Pf(x'(0)-h))+常数
%output
%Kf:控制增益
%g:控制常量
%Pf: 代价矩阵
%h: 目标函数偏移项
%input
%A：状态转移矩阵
%B: 控制矩阵
%f: 常数
%Q，R，F: 代价矩阵
%N:控制步长
%x0:初始状态
%z:追踪定点
    if(rank(ctrb(A,B)) ~= size(x0,1))
        error('(A,B) 必须可控');
    end
    K = zeros(size(B,2),size(A,1),N);
    P = zeros(size(F,1),size(F,2),N+1);
    v = zeros(size(A,1),1,N+1);
    P(:,:,N+1) = F;
    v(:,:,N+1) = z - f;
    for i=N:-1:1
        K(:,:,i) = (B'*P(:,:,i+1)*B+R)\B'*P(:,:,i+1)*A;
        P(:,:,i) = (A-B*K(:,:,i))'*P(:,:,i+1)*(A-B*K(:,:,i))+...
            K(:,:,i)'*R*K(:,:,i)+Q;
        v(:,:,i) = inv(P(:,:,i))*(((A-B*K(:,:,i))'*P(:,:,i+1)*(A-B*K(:,:,i))+...
            K(:,:,i)'*R*K(:,:,i))*inv(A)+Q)*v(:,:,i+1);
    end
    Kf = -K(:,:,1);
    g = -A\v(:,:,2);
    Pf = P(:,:,1);
    h = v(:,:,1);
end
%end dlqr_N