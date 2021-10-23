function [Xf,Pf,h] = TerminalSet_and_Penalty(sys,z)
%计算系统sys的关于定点z的不变集
%sys-------x(k+1) = Ax(k)+Bu(k)+f
%使用mpc控制
%----------u = Kf(x+d)----------
%计算  S(k+1) = {x|(A+B*Kf)x+B*Kf*d+f属于S(k)}
[Kf,d,Pf,h] = mpc_Gain_Penalty(sys,z);
model2 = LTISystem('A',sys.A+sys.B*Kf,'f',sys.B*Kf*d+sys.f);
model2.x.min = sys.x.min;
model2.x.max = sys.x.max;
model2.x.penalty = sys.x.penalty;
Xf = model2.invariantSet();
U0 = Polyhedron('lb',sys.u.min,'ub',sys.u.max);
if(size(Kf,1) == size(Kf,2))
    cond2 = U0.invAffineMap(Kf,Kf*d);
else
    cond2=Polyhedron('A',[-Kf;Kf],'b',[U0.b(2,1)-Kf*d;-(-U0.b(1,1)-Kf*d)]);%注意修改控制约束
end
Xf = Xf.intersect(cond2).minHRep();
end