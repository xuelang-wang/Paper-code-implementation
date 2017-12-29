%On clock synchronization algorithm for wireless sensor network under unknowm delay
%autor:Mei Leng and Yik-Chung Wu
%2016.12.23 by xu
%==========================================================================
close all;
%==========================================================================
%%
clear;
clc;
N = 11;%measures times
h = 1;
for SNR = 0:15:30
    %parameters setting
    O = 30;%simulation runs times
    H = 25;%中心节点发送时间间隔
    w_sigma = 0.3 * H;%中心节点发送时刻噪声的方差
    G = 30;%参考节点发送时间间隔
    v_sigma = 0.3 * G;%参考节点发送时刻噪声的方差
    d_lb = 0 + eps;%lower bound of fixed delay(open interval)
    d_ub = 10;%upper bound of fixed delay
    beta1_lb = 0.9;%lower bound of relative skew
    beta1_ub = 1.1;%upper bound of relative skew
    beta0_lb = -10;%lower bound of relative offset
    beta0_ub = 10;%upper bound of  relative offset
    sigma = (H^2+G^2)/(power(10,SNR/10));%variance of random delay Xi and Yi
    k = 1;
    for alpha = 1:10%measures times
        %time stamps
        T_1 = zeros(N,O);
        T_2 = zeros(N,O);
        T_3 = zeros(N,O);
        T_4 = zeros(N,O);
        d = zeros(1,O);%real fixed delay
        beta1 = zeros(1,O);%real clock skew
        beta0 = zeros(1,O);%real clock offset
        beta1_mle = zeros(1,O);%estimates clock skew
        beta0_mle = zeros(1,O);%estimates clock offset
        for o = 1:O
            d(1,o) = unifrnd(d_lb,d_ub);
            beta1(1,o) = unifrnd(beta1_lb,beta1_ub);
            beta0(1,o) = unifrnd(beta0_lb,beta0_ub);
            %produce measures
            for i = 1:N
                T_1(i,o) = i * H + normrnd(0,sqrt(w_sigma));
                T_3(i,o) = i * G + normrnd(0,sqrt(v_sigma));
                T_2(i,o) = beta1(1,o) * T_1(i,o) + beta0(1,o) + beta1(1,o) * (d(1,o) + normrnd(0,sqrt(sigma)));
                T_4(i,o) = (T_3(i,o) - beta0(1,o)) / beta1(1,o) + d(1,o) + normrnd(0,sqrt(sigma));
            end
            %MLLE <2>alpha= alpha star 
            %produce alpha star
            j = mod(N,3);
            n = (N - j)/3;
            alpha_star_star = 2*n+ceil(j/2);
            for alpha_star = 1:10
                for m = 1:N-alpha_star
                    D_1(m,alpha_star,o) = T_1(alpha_star+m,o)-T_1(m,o);
                    D_2(m,alpha_star,o) = T_2(alpha_star+m,o)-T_2(m,o);
                    D_3(m,alpha_star,o) = T_3(alpha_star+m,o)-T_3(m,o);
                    D_4(m,alpha_star,o) = T_4(alpha_star+m,o)-T_4(m,o);
                end
                beta1_mlle_PB2(alpha_star,o) = 2*sigma*beta1(1,o)^4/...
                    (sum(beta1(1,o)^2*D_1(:,alpha_star,o).^2+D_3(:,alpha_star,o).^2+6*beta1(1,o)^2*sigma));
            end
        end
        %MLLE  <2>alpha = alpha star
        beta1_mlle_PB2_mes(h,k) = sum(beta1_mlle_PB2(k,:),2)/O;
        k = k + 1;
    end
    figure(2)
    semilogy(1:10,beta1_mlle_PB2_mes(h,:),'ok-');
    hold on;
    semilogy(alpha_star_star,beta1_mlle_PB2_mes(h,alpha_star_star),'rx','MarkerSize',20);
    gtext(['SNR = ',num2str(SNR),' bB']);
    h=h+1;
end
figure(2)
grid on;
legend('N = 11');
xlabel('\alpha');
ylabel('PB_g of the clock skew(\beta_1)');
%==========================================================================
%%
clear;
clc;
%parameters setting
O = 10000;%simulation runs times
H = 25;%中心节点发送时间间隔
w_sigma = 0.3 * H;%中心节点发送时刻噪声的方差
G = 30;%参考节点发送时间间隔
v_sigma = 0.3 * G;%参考节点发送时刻噪声的方差
d_lb = 0 + eps;%lower bound of fixed delay(open interval)
d_ub = 10;%upper bound of fixed delay
beta1_lb = 0.9;%lower bound of relative skew
beta1_ub = 1.1;%upper bound of relative skew
beta0_lb = -10;%lower bound of relative offset
beta0_ub = 10;%upper bound of  relative offset
SNR = 30;%dB
sigma = (H^2+G^2)/(power(10,SNR/10));%variance of random delay Xi and Yi

k = 1;
for N = 6:2:18%measures times
    %time stamps
    T_1 = zeros(N,O);
    T_2 = zeros(N,O);
    T_3 = zeros(N,O);
    T_4 = zeros(N,O);
    d = zeros(1,O);%real fixed delay
    beta1 = zeros(1,O);%real clock skew
    beta0 = zeros(1,O);%real clock offset
    beta1_mle = zeros(1,O);%estimates clock skew
    beta0_mle = zeros(1,O);%estimates clock offset
    beta1_CRLB = zeros(1,O);
    beta0_CRLB = zeros(1,O);
    d_CRLB = zeros(1,O);
    beta1_mlle_1 = zeros(1,O);
    beta0_mlle_1 = zeros(1,O);
    beta1_mlle_PB1 = zeros(1,O);
    beta0_mlle_PB1 = zeros(1,O);
    beta1_mlle_2 = zeros(1,O);
    beta0_mlle_2 = zeros(1,O);
    beta1_mlle_PB2 = zeros(1,O);
    beta0_mlle_PB2 = zeros(1,O);
    beta1_my = zeros(1,O);
    beta0_my = zeros(1,O);
    beta1_my_PB = zeros(1,O);
    beta0_my_PB = zeros(1,O);
    for o = 1:O
        d(1,o) = unifrnd(d_lb,d_ub);
        beta1(1,o) = unifrnd(beta1_lb,beta1_ub);
        beta0(1,o) = unifrnd(beta0_lb,beta0_ub);
        for i = 1:N
            T_1(i,o) = i * H + normrnd(0,sqrt(w_sigma));
            T_3(i,o) = i * G + normrnd(0,sqrt(v_sigma));
            T_2(i,o) = beta1(1,o) * T_1(i,o) + beta0(1,o) + beta1(1,o) * (d(1,o) + normrnd(0,sqrt(sigma)));
            T_4(i,o) = (T_3(i,o) - beta0(1,o)) / beta1(1,o) + d(1,o) + normrnd(0,sqrt(sigma));
        end

        %MLE
        %matrix form
        Ts = [T_1(:,o);-T_4(:,o)];
        Tp = [T_2(:,o) -ones(N,1);-T_3(:,o) ones(N,1)];
        P = eye(2*N) - Tp * inv(Tp' * Tp) * Tp';
        d_est = -(ones(1,2*N) * P * Ts + Ts' * P * ones(2*N,1))/(2 * ones(1,2*N) * P * ones(2*N,1));
        theta_est = inv(Tp' * Tp)*Tp'*(Ts + d_est * ones(2*N,1));
        beta1_mle(1,o) = 1/theta_est(1,1);
        beta0_mle(1,o) = theta_est(2,1)/theta_est(1,1);
        %CRLB
        A = beta1(1,o)*beta1(1,o)*(T_1(:,o)+d(1,o))'*(T_1(:,o)+d(1,o)) +...
            N*beta1(1,o)^2*sigma + (T_3(:,o)-beta0(1,o))'*(T_3(:,o)-beta0(1,o))/(beta1(1,o)^4);
        B = sum(beta1(1,o)*(T_1(:,o)+d(1,o)) + (T_3(:,o)-beta0(1,o)))/(beta1(1,o)^3);
        C = sum(beta1(1,o)*(T_1(:,o)+d(1,o))-(T_3(:,o)-beta0(1,o)))/(beta1(1,o)^2);
        beta1_CRLB(1,o) = 2 * N *sigma/(2*N*A-beta1(1,o)^2*B^2-C^2);
        beta0_CRLB(1,o) = sigma*beta1(1,o)^2*(2*N*A-C^2)/(2*N*(2*N*A-beta1(1,o)^2*B^2-C^2));
        d_CRLB(1,o) = sigma*(2*N*A-beta1(1,o)^2*B^2)/(2*N*(2*N*A-beta1(1,o)^2*B^2-C^2));
        %MLLE <1>alpha = N-1
        alpha = N - 1;
        W_1 = T_1(N,o) - T_1(1,o);
        W_2 = T_2(N,o) - T_2(1,o);
        W_3 = T_3(N,o) - T_3(1,o);
        W_4 = T_4(N,o) - T_4(1,o);
        beta1_mlle_1(1,o) = (W_2^2+W_3^2)/(W_1*W_2+W_4*W_3);
        beta0_mlle_1(1,o) = sum(T_2(:,o) + T_3(:,o)-beta1_mlle_1(1,o)*(T_1(:,o)+T_4(:,o)))/(2*N);
        beta1_mlle_PB1(1,o) = 2*sigma*beta1(1,o)^4/(beta1(1,o)^2*W_1^2+W_3^2+6*beta1(1,o)^2*sigma);
        beta0_mlle_PB1(1,o) = sigma*beta1(1,o)^2/(2*N)+beta1_mlle_PB1(1,o)/...
            (4*N^2)*(sum(beta1(1,o)*(T_1(:,o)+d(1,o))+(T_3(:,o)-beta0(1,o)))*...
            sum(beta1(1,o)*(T_1(:,o)+d(1,o))+(T_3(:,o)-beta0(1,o)))/(beta1(1,o)^2)+N*sigma);
        
        %MLLE <2>alpha= alpha star 
        %produce alpha star
        j = mod(N,3);
        n = (N - j)/3;
        alpha_star = 2*n+ceil(j/2);
        for m = 1:N-alpha_star
            D_1(m,o) = T_1(alpha_star+m,o)-T_1(m,o);
            D_2(m,o) = T_2(alpha_star+m,o)-T_2(m,o);
            D_3(m,o) = T_3(alpha_star+m,o)-T_3(m,o);
            D_4(m,o) = T_4(alpha_star+m,o)-T_4(m,o);
        end
        beta1_mlle_2(1,o) = (D_2(:,o)'*D_2(:,o)+D_3(:,o)'*D_3(:,o))/(D_1(:,o)'*D_2(:,o)+D_4(:,o)'*D_3(:,o));
        beta0_mlle_2(1,o) = sum(T_2(:,o) + T_3(:,o)-beta1_mlle_2(1,o)*(T_1(:,o)+T_4(:,o)))/(2*N);
        beta1_mlle_PB2(1,o) = 2*sigma*beta1(1,o)^4/(sum(beta1(1,o)^2*D_1(:,o).^2+D_3(:,o).^2+6*beta1(1,o)^2*sigma));
        beta0_mlle_PB2(1,o) = sigma*beta1(1,o)^2/(2*N)+beta1_mlle_PB2(1,o)/...
            (4*N^2)*(sum(beta1(1,o)*(T_1(:,o)+d(1,o))+(T_3(:,o)-beta0(1,o)))*sum(beta1(1,o)*(T_1(:,o)+d(1,o))+(T_3(:,o)-beta0(1,o)))/(beta1(1,o)^2)+N*sigma);
        
        %Proposed algorithm
        Ts_1 = T_1(:,o)+T_4(:,o);
        Tp_1 = [T_2(:,o)+T_3(:,o) -2*ones(N,1)];
        theta_1 = inv(Tp_1'*Tp_1)*Tp_1'*Ts_1;
        beta1_my(1,o) = 1/theta_1(1,1);
        beta0_my(1,o) = theta_1(2,1)/theta_1(1,1);
        K = sum((T_1(:,o) + d(1,o) + (T_3(:,o)-beta0(1,o)/beta1(1,o))).^2+3*sigma)/beta1(1,o)^2;
        beta1_my_PB(1,o) = 2*N*sigma/(N*K-beta1(1,o)^2*B^2);
        beta0_my_PB(1,o) = sigma*beta1(1,o)^2*K/(2*N*K-2*beta1(1,o)^2*B^2);
    end
    %MLE
    beta1_mle_mes(1,k) = (beta1 - beta1_mle)*(beta1 - beta1_mle)'/O;
    beta0_mle_mes(1,k) = (beta0 - beta0_mle)*(beta0 - beta0_mle)'/O;
    beta1_CRLB_mes(1,k) = sum(beta1_CRLB,2) / O;
    beta0_CRLB_mes(1,k) = sum(beta0_CRLB,2) / O;
    %MLLE  <1>alpha = N-1
    beta1_mlle_mes1(1,k) = (beta1 - beta1_mlle_1)*(beta1 - beta1_mlle_1)'/O;
    beta0_mlle_mes1(1,k) = (beta0 - beta0_mlle_1)*(beta0 - beta0_mlle_1)'/O;
    beta1_mlle_PB1_mes(1,k) = sum(beta1_mlle_PB1,2)/O;
    beta0_mlle_PB1_mes(1,k) = sum(beta0_mlle_PB1,2)/O;
    %MLLE  <2>alpha = alpha star
    beta1_mlle_mes2(1,k) = (beta1 - beta1_mlle_2)*(beta1 - beta1_mlle_2)'/O;
    beta0_mlle_mes2(1,k) = (beta0 - beta0_mlle_2)*(beta0 - beta0_mlle_2)'/O;
    beta1_mlle_PB2_mes(1,k) = sum(beta1_mlle_PB2)/O;
    beta0_mlle_PB2_mes(1,k) = sum(beta0_mlle_PB2)/O;
    %proposed algorithm
    beta1_my_mes(1,k) = (beta1 - beta1_my)*(beta1 - beta1_my)'/O;
    beta0_my_mes(1,k) = (beta0 - beta0_my)*(beta0 - beta0_my)'/O;
    beta1_my_PB_mes(1,k) = sum(beta1_my_PB)/O;
    beta0_my_PB_mes(1,k) = sum(beta0_my_PB)/O;
    
    k = k + 1;
end
t = 6:2:18;
figure(3);
hold on;

plot(t,beta0_mlle_mes1,'sk-',t,beta0_mlle_PB1_mes,'dr--');
plot(t,beta0_mlle_mes2,'dk-',t,beta0_mlle_PB2_mes,'+r--');
plot(t,beta0_my_mes,'<k-',t,beta0_my_PB_mes,'xr--');
plot(t,beta0_mle_mes,'>k-',t,beta0_CRLB_mes,'or--');
legend('GE-(\alpha=N-1)','PB_g-(\alpha=N-1)','GE-(\alpha=\alpha^*)',...
    'PB_g-(\alpha=\alpha^*)','Proposed','PB_p','MLE','CRLB');
grid on;
xlabel('N');
ylabel('MSE of the clock offset(\beta_0)');
axis([6,18,0.2,1]);
figure(4);

semilogy(t,beta1_mlle_mes1,'sk-',t,beta1_mlle_PB1_mes,'dr--');
hold on;
semilogy(t,beta1_mlle_mes2,'dk-',t,beta1_mlle_PB2_mes,'+r--');
plot(t,beta1_my_mes,'<k-',t,beta1_my_PB_mes,'xr--');
semilogy(t,beta1_mle_mes,'>k-',t,beta1_CRLB_mes,'or--');
legend('GE-(\alpha=N-1)','PB_g-(\alpha=N-1)','GE-(\alpha=\alpha^*)',...
    'PB_g-(\alpha=\alpha^*)','Proposed','PB_p','MLE','CRLB','location','Best');
grid on;
xlabel('N');
ylabel('MSE of the clock skew (\beta_1)')
axis([6,18,1e-6,1e-4]);


