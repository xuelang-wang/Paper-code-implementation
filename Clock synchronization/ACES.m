%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ACES_Adaptive clock estimation and synchronization using kalman filtering
%author: Benjamin R.Hamiton
%years: 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
%%
%---------------------Figure3:Clock Skew Tracking Example------------------
%Skew estimation with a scalar kalman filter
clear
clc
%parameter settings
tau_0 = 222;%sample period
v = 3600;%the normlizer to model different decaying rates
rho = 0.998;%position number close to 1
p = power(rho,tau_0/v);%AR(1)parameter
sigma_alpha = 0.01;%standard deviation of skew???
sigma_eta(1,1) = sqrt((1 - p^2))*sigma_alpha;%standard of process;
sigma_eta(1,2) = 0.1*sqrt((1 - p^2))*sigma_alpha;%standard deviation of process;
sigma_eta(1,3) = 10*sqrt((1 - p^2))*sigma_alpha;%standard deviation of process;
% sigma_v = sqrt(2)*sigma_w;%standard deviation of measurement noise;
sigma_v = 0.05;%affect speed of converges???
N = 50;
O = 20;
alpha_real = zeros(N,3);
alpha_est = zeros(N,3);
M = zeros(N,3);
%intialization(need to adjust)
for o = 1:O
    alpha_real(1,:,o) = normrnd(0,sigma_alpha,1,3);%????
    alpha_est(1,:,o) = alpha_real(1,:,o)+normrnd(0,sigma_v,1,3)/tau_0;%???
    M(1,1) = sigma_alpha^2;
    M(1,2) = sigma_alpha^2;
    M(1,3) = sigma_alpha^2;
    for i = 1:N-1
        alpha_real(i+1,:,o) = p*alpha_real(i,:,o)+normrnd(0,sigma_eta(1,:));%real update
        deta_theta = alpha_real(i+1,:,o)*tau_0+normrnd(0,sigma_v,1,3);%produce measurement
        temp = p^2*M(i,:) + sigma_eta(1,:).^2;%predictive MSE
        Gain = tau_0*temp./(sigma_v^2+tau_0^2.*temp);%kalman gain
        alpha_est(i+1,:,o) = p*alpha_est(i,:,o) + Gain.*(deta_theta - tau_0*p*alpha_est(i,:,o));%estimation update
        M(i+1,:) = (1-tau_0*Gain).*temp;%MMSE
    end
end
alpha_real_ave = zeros(N,3);
alpha_est_ave = zeros(N,3);
for o = 1:O
    alpha_real_ave(:,:) = alpha_real_ave(:,:) + alpha_real(:,:,o);
    alpha_est_ave(:,:) = alpha_est_ave(:,:) + alpha_est(:,:,o);
end
alpha_real_ave = alpha_real_ave/O;
alpha_est_ave = alpha_est_ave/O;
figure
n = 100:222:222*N;
plot(n,alpha_real_ave(:,1),'bx-');
hold on;
plot(n,alpha_est_ave(:,1),'g+-',n,alpha_est_ave(:,2),'ro-',n,alpha_est_ave(:,3),'c.-');
legend('\alpha','Estimated\alpha(\sigma_{\eta})','Estimated\alpha(0.1\sigma_{\eta})','Estimated\alpha(10\sigma_{\eta})')
figure
semilogy(n,M(:,1),'b+-',n,M(:,2),'go-',n,M(:,3),'r.-');
legend('Estimated\alpha(\sigma_{\eta})','Estimated\alpha(0.1\sigma_{\eta})','Estimated\alpha(10\sigma_{\eta})')

%%
%---------------------------figure4:(a)------------------------------------ 
clear
clc
%parameter settings
N = 41;
tau_0 = zeros(N,1);
for i = 1:N
    tau_0(i,1) = 1e-5*power(10,0.5*(i-1));
end
sigma_alpha_1 = 0.0001;
sigma_inf_1=zeros(N,6);
rho_1 = 1-1e-11;
v_1 = 3600;
sigma_v_1 = [0 0.1 1 10 100 1000];
for i = 1:length(sigma_v_1)
    sigma_inf_1(:,i) = (1- power(rho_1,2*tau_0/v_1)).*(sigma_alpha_1^2*tau_0.^2-sigma_v_1(1,i)^2)./(2*tau_0.^2)+...
        sqrt((1- power(rho_1,2*tau_0/v_1)).^2.*(sigma_alpha_1^2*tau_0.^2-sigma_v_1(1,i)^2).^2+...
        4*(1- power(rho_1,2*tau_0/v_1)).*tau_0.^2*sigma_v_1(1,i)^2*sigma_alpha_1^2)./(2*tau_0.^2);
end

figure
loglog(tau_0,sigma_inf_1(:,1),'b-',tau_0,sigma_inf_1(:,2),'g-x',tau_0,sigma_inf_1(:,3),'r-+',...
    tau_0,sigma_inf_1(:,4),'-c.',tau_0,sigma_inf_1(:,5),'m-s',tau_0,sigma_inf_1(:,6),'y-d');
legend('\sigma_v=0','\sigma_v=0.1','\sigma_v=1','\sigma_v=10','\sigma_v=100','\sigma_v=1000','location','best');
xlabel('\tau_0(s)');
ylabel('\Sigma_{\infty}([s/s]^2)');
axis([1e-5,1e15,1e-24,1e-8]);

%---------------------------figure4:(b)------------------------------------ 
%parameter settings
sigma_alpha_2 = [0,1e-6,1e-5,1e-4,1e-3,0.01];
sigma_inf_2=zeros(N,6);
rho_2 = 1-1e-11;
v_2 = 3600;
sigma_v_2 = 10;
for i = 1:length(sigma_alpha_2)
    sigma_inf_2(:,i) = (1- power(rho_2,2*tau_0/v_2)).*(sigma_alpha_2(1,i)^2*tau_0.^2-sigma_v_2^2)./(2*tau_0.^2)+...
        sqrt((1- power(rho_2,2*tau_0/v_2)).^2.*(sigma_alpha_2(1,i)^2*tau_0.^2-sigma_v_2^2).^2+...
        4*(1- power(rho_2,2*tau_0/v_2)).*tau_0.^2*sigma_v_2^2*sigma_alpha_2(1,i)^2)./(2*tau_0.^2);
end

figure
loglog(tau_0,sigma_inf_2(:,1),'b-',tau_0,sigma_inf_2(:,2),'g-x',tau_0,sigma_inf_2(:,3),'r-+',...
    tau_0,sigma_inf_2(:,4),'-c.',tau_0,sigma_inf_2(:,5),'m-s',tau_0,sigma_inf_2(:,6),'y-d');
legend('\sigma_{\alpha}=0','\sigma_{\alpha}=1e-006','\sigma_{\alpha}=1e-005',...
    '\sigma_{\alpha}=0.0001','\sigma_{\alpha}=0.001','\sigma_{\alpha}=0.01','location','best');
xlabel('\tau_0(s)');
ylabel('\Sigma_{\infty}([s/s]^2)');
axis([1e-5,1e15,1e-18,1e-4]);
    
%---------------------------figure4:(c)------------------------------------   
%parameter settings
sigma_alpha_3 = 0.0001;
sigma_inf_3=zeros(N,6);
rho_3 = [1,1-1e-15,1-1e-13,1-1e-11,1-1e-9,1-1e-7];
v_3 = 3600;
sigma_v_3 = 10;
for i = 1:length(rho_3)
    sigma_inf_3(:,i) = (1- power(rho_3(1,i),2*tau_0/v_3)).*(sigma_alpha_3^2*tau_0.^2-sigma_v_3^2)./(2*tau_0.^2)+...
        sqrt((1- power(rho_3(1,i),2*tau_0/v_3)).^2.*(sigma_alpha_3^2*tau_0.^2-sigma_v_3^2).^2+...
        4*(1- power(rho_3(1,i),2*tau_0/v_3)).*tau_0.^2*sigma_v_3^2*sigma_alpha_3^2)./(2*tau_0.^2);
end

figure
loglog(tau_0,sigma_inf_3(:,1),'b-',tau_0,sigma_inf_3(:,2),'g-x',tau_0,sigma_inf_3(:,3),'r-+',...
    tau_0,sigma_inf_3(:,4),'-c.',tau_0,sigma_inf_3(:,5),'m-s',tau_0,sigma_inf_3(:,6),'y-d');
legend('\rho=1-0','\rho=1-1e-015','\rho=1-1e-013','\rho=1-1e-011','\rho=1-1e-009','\rho=1-1e-007','location','best');
xlabel('\tau_0(s)');
ylabel('\Sigma_{\infty}([s/s]^2)');
axis([1e-5,1e15,1e-18,1e-8]);

%%    
%-----------figure5:simulated Sigma_inf for Skew Tracking------------------     
clear
clc
%parameter settings
O = 50;%synchronization period
M = 50;%runs times,,average no effect
N = 41;
tau_0 = zeros(N,1);
for i = 1:N
    tau_0(i,1) = 1e-5*power(10,0.5*(i-1));
end
sigma_alpha = 0.0001;
sigma_inf=ones(N,6);
rho = 1-1e-11;
v = 3600;
sigma_v = [0 0.1 1 10 100 1000];
for i = 1:length(sigma_v)
    for  o = 1:O
        sigma_inf(:,i) = power(rho,2*tau_0/v).*sigma_v(1,i)^2.*sigma_inf(:,i)...
            ./(sigma_v(1,i)^2+tau_0.^2.*sigma_inf(:,i))+(1- power(rho,2*tau_0/v))*sigma_alpha^2;
    end
end

figure
loglog(tau_0,sigma_inf(:,1),'b-',tau_0,sigma_inf(:,2),'g-x',tau_0,sigma_inf(:,3),'r-+',...
    tau_0,sigma_inf(:,4),'-c.',tau_0,sigma_inf(:,5),'m-s',tau_0,sigma_inf(:,6),'y-d');
legend('\sigma_v=0','\sigma_v=0.1','\sigma_v=1','\sigma_v=10','\sigma_v=100','\sigma_v=1000','location','best');
xlabel('\tau_0(s)');
ylabel('\Sigma_{\infty}([s/s]^2)');
axis([1e-5,1e15,1e-24,1e-8]);
%%
%Skew and offset estimation with a vector kalman filter
clear
clc
%parameter settings
tau_0 = 1;%sample period
v = 3600;%the normlizer to model different decaying rates
rho = 1 - 2e-6;%position number close to 1
p = power(rho,tau_0/v);%AR(1)parameter
sigma_alpha = 0.0001;%standard deviation of skew???
sigma_eta = sqrt((1 - p^2))*sigma_alpha;%standard of process;
sigma_v = sqrt(1e-6);%affect speed of converges
A = [1 tau_0;0 p];
b = [1 0]';
N = 30;
O = 20;
%intialization(need to adjust)
x_real = zeros(2,N+1);
x_est = zeros(2,N+1);
theta_ob = zeros(1,N+1);%observed value
theta_alpha = zeros(1,N+1);%alpha-estimated
M = zeros(2,2,N+1);
M(:,:,1) = diag([sigma_v^2,sigma_alpha^2]);
alpha_real = zeros(1,N+1);
alpha_est = zeros(1,N+1);

M_alpha = zeros(N+1,1);
M_alpha(1,1) = sigma_alpha^2;
M_theta = zeros(N+1,1);
M_theta(1,1) = sigma_v^2;
for i = 1:N
    x_real(:,i+1) = A*x_real(:,i)+[0;normrnd(0,sigma_eta)];%real update
    theta_ob(1,i+1) = b'*x_real(:,i+1)+normrnd(0,sigma_v);%produce measurement
    MSE = A*M(:,:,i)*A'+diag([0,sigma_eta^2]);%predictive MSE
    Kalman_Gain = MSE*b*inv(sigma_v^2+b'*MSE*b);%kalman gain
    x_est(:,i+1) = A*x_est(:,i)+Kalman_Gain*(theta_ob(1,i+1)-b'*A*x_est(:,i)); 
    M(:,:,i+1) = (eye(2,2) - Kalman_Gain*b')*MSE;%MMSE;
    
    alpha_real(1,i+1) = p*alpha_real(1,i)+normrnd(0,sigma_eta);%real update
    deta_theta = alpha_real(1,i+1)*tau_0+normrnd(0,sigma_v);%produce measurement
    temp = p^2*M_alpha(i,:) + sigma_eta^2;%predictive MSE
    Gain = tau_0*temp/(sigma_v^2+tau_0^2*temp);%kalman gain
    alpha_est(1,i+1) = p*alpha_est(1,i) + Gain*(deta_theta - tau_0*p*alpha_est(1,i));%estimation update
    theta_alpha(1,i+1) = theta_alpha(1,i) + alpha_est(1,i)*tau_0;
    M_alpha(i+1,:) = (1-tau_0*Gain)*temp;%MMSE
    M_theta(i+1,:) = M_theta(i,:) + tau_0^2*M_alpha(i,:)+sigma_v^2;
end
for i = 1:N+1
    M_est(1,i) = M(1,1,i);
    M_ob(1,i) = (theta_ob(1,i)-x_real(1,i))^2;
end
figure
n = 0:N;
plot(n,x_real(1,:),'bx-',n,x_est(1,:),'g-o',n,theta_alpha,'-r+',n,theta_ob,'c*-');
axis([0,30,-8e-3,6e-3]);
legend('True\theta[n]','Estimated\theta[n]','\alpha-estimated\theta[n]','Observed\theta[n]');
figure
semilogy(n,M_est(1,:),'go-',n(1,2:end),M_theta(2:end,1),'r-+',n,M_ob(1,:),'c-*');
legend('Estimated\theta[n]','\alpha-estimated\theta[n]','Observed\theta[n]','location','best');

%%
%------------Figure7:Simulated MSEs vs. time-------------------------------
clear
clc
%parameter settings
tau_0 = zeros(11,1);%sample period
for i = 1:11
    tau_0(i,1) = 3*power(10,i-6);
end
v = 3600;%the normlizer to model different decaying rates
rho = 1 - 2e-6;%position number close to 1
p = power(rho,tau_0/v);%AR(1)parameter
sigma_alpha = 0.0002;%standard deviation of skew???
sigma_eta = sqrt((1 - p.^2))*sigma_alpha;%standard of process;
sigma_v = sqrt(1e-5);%affect speed of converges
for i = 1:11
    A = [1 tau_0(i,1);0 p(i,1)];
end
b = [1 0]';
N = 10000;%100000 synchronization period
for j = 1:11
    MSE(:,:,1) = [1e-6 0;0 1e-8];%intialization(need to adjust)
    MSE_skew(j,1) = MSE(2,2,1);
    MSE_offset(j,1) = MSE(1,1,1);
    for i = 1:N
        Kalman_Gain = MSE(:,:,i)*b*inv(sigma_v^2+b'*MSE(:,:,i)*b);%kalman gain
        M = (eye(2,2) - Kalman_Gain*b')*MSE(:,:,i);%MMSE;
        MSE(:,:,i+1) = A*M*A'+diag([0,sigma_eta(j,1)^2]);%predictive MSE
        MSE_skew(j,i+1) = MSE(2,2,i+1);
        MSE_offset(j,i+1) = MSE(1,1,i+1);
    end
    t = power(10,-4.5):tau_0(j,1):(N+5)*tau_0(j,1);
    figure(9)
    loglog(t(1:N+1),MSE_skew(j,:));
    hold on;
%     axis([1e-5,1e15,1e-18,1e-6]);
    figure(10)
    loglog(t(1:N+1),MSE_offset(j,:));
    hold on;
%     axis([1e-5,1e15,1e-10,1e-5]);
end
% figure(9)
% legend('True\theta[n]','Estimated\theta[n]','\alpha-estimated\theta[n]','Observed\theta[n]');
% figure(10)
% legend('Estimated\theta[n]','\alpha-estimated\theta[n]','Observed\theta[n]','location','best');














