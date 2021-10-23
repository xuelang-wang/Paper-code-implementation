function [S,flag] = invariant_error(AK,W)
%system
%------e(k+1) = AK*e(k)+w(k)--------
%--------w(k)属于紧集W--------------
%-----AK的特征值全部小于1，则存在一个集合S
%-----AK*S+W属于S,'+'为闵可夫斯基加
flag = 0;
iter = 50;
S = W;
for i = 1:iter
    S1 = AK*S+W;
    S1.minHRep();
    if(S.eq(S1))
        flag = 1;
        break;
    end
    S = S1;
end
if(~flag)
    warning('Computation finished without convergence.');
end
