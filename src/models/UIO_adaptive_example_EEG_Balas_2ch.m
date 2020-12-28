%% Unknown Input Observer with Adaptive State Estimation
% Implementation from the following paper:
% B. Alenezi, J. Hu and S. H. {ak, "Adaptive Unknown Input and State Observers," 2019 American Control Conference (ACC), Philadelphia, PA, USA, 2019, pp. 2434-2439, doi: 10.23919/ACC.2019.8815288.
%% Setup
%opengl software
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%% Prepare Data
fs=512;
dt=inv(fs);
EEG=load('/home/tris/Github/ardent_adaptive_UIO/data/processed/group_01_sujet_01.mat');
%Y_data=EEG.samples(:,[2,8,16,20,24,30]);
Y_data=EEG.samples(:,[2,30]);
t=EEG.samples(:,1);
Y_data=detrend(Y_data);
Y_data=bandpass(Y_data,[4 45],fs);
Y_data=Y_data';

    
Y_sim=array2timetable(Y_data(:,1:512*60)','SampleRate',fs);
Y_sim1=Y_sim(:,1);
Y_sim2=Y_sim(:,2);
% Y_sim3=Y_sim(:,3);
% Y_sim4=Y_sim(:,4);
% Y_sim5=Y_sim(:,5);
% Y_sim6=Y_sim(:,6);
%% OMA
order = 10;
s = 2*order;
opt_order=4;
[A_data,C_data,G_data,R0_data] = ssidata(Y_data,order,s);
err = [0.01,0.05,0.98];
[IDs_cov] = plotstab(A_data,C_data,Y_data,dt,[],err);
%% Vectors and Analysis
[fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,dt);
Phi_data=Phi_data';
norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest
ss_d=ss(A_data{opt_order},[],C_data{opt_order},[],fs);
ss_c=d2c(ss_d);

B2=ones(opt_order,1);


%% Wronskian test
syms t
f = [1 cos(t) sin(t)];
W = [f; f; f];
for n=2:3
    W(n,:) = diff(W(n-1,:)); 
end
f=f.'
%% Double Integrator with Unknown Input
syms t
epsilon=4/5;
Lx=5;
Am=[0 1;-5 -4];
A=[0 1;-5+Lx -4+epsilon*Lx]
B=[0;1];
C=[1, epsilon];
Phi_k = [1 cos(t) sin(t)].';
gamma_e=1;
gamma_u=2;
g1=4;
g2=-0.3;
g3=0.1;
ThetaLstar=[g1,g2,g3];
%% Sim
sim_out=sim('UIO_balas_01.slx',60);
figure
ax1=subplot(2,3,1);
plot(sim_out.y)
hold on
plot(sim_out.yhat)
grid on
%xlim([18,22])
legend('$y$','$\hat{y}$')
title('($y$ and $\hat{y}$)')
xlabel('Time (s)')
ylabel('Output State ($y$)')
ax2=subplot(2,3,2);
plot(sim_out.ey)
grid on
legend('$e_y$')
title('State Error')
xlabel('Time (s)')
ylabel('Error ($e_y$)')

ax3=subplot(2,3,3);
plot(sim_out.u)
hold on
plot(sim_out.uhat)
grid on
title('Input State ($u$ and $\hat{u}$)')
xlabel('Time (s)')
ylabel('Input State')
legend('$u$','$\hat{u}$')
ax4=subplot(2,3,4);
plot(sim_out.eu)
grid on

title('Input Error ($e_u$)')
xlabel('Time (s)')
ylabel('Input Error ($e_u$)')
legend('$e_u$')

ax4=subplot(2,3,5);
plot(sim_out.ThetaL)
grid on

title('Basis $\alpha$')
xlabel('Time (s)')
ylabel('$\alpha$')
legend('$\alpha_{1}$','$\alpha_{\sin}$','$\alpha_{\cos}$')

%linkaxes([ax1,ax2],'y');
sgtitle(['UIO for Input of ', num2str(g2),'$\sin{t}+$',num2str(g3),'$\cos{t}+$',num2str(g1)])


