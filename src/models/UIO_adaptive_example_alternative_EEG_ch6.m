%% Unknown Input Observer with Adaptive State Estimation
% Implementation from our own algorithm on EEG
%% Setup
%opengl software
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%% Prepare Data
fs=512;
dt=inv(fs);
EEG=load('F:\Documents\MATLAB\group_01_mat\group_01_sujet_01.mat');
Y_data=EEG.samples([100:end],[2,8,16,20,24,30]);
%Y_data=EEG.samples([100:end],[2]);

t=EEG.samples([100:end],1);
Y_data=detrend(Y_data);
Y_data=bandpass(Y_data,[1 49],fs); % PLAY with this
Y_data=Y_data';

    
Y_sim=array2timetable(Y_data(:,1:512*60*3)','SampleRate',fs);
Y_sim1=Y_sim(:,1);
Y_sim2=Y_sim(:,2);
Y_sim3=Y_sim(:,3);
Y_sim4=Y_sim(:,4);
Y_sim5=Y_sim(:,5);
Y_sim6=Y_sim(:,6);
%% OMA
order = 40;
s = 2*order;
opt_order=6;
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

B=ones(opt_order,1);

%%

Am=ss_c.A;
C=ss_c.C;
mats=[];
for i=1:1:40
    mat=[0 1;-i^2 0]
    mats=blkdiag(mats,mat)
end
%Fu=[0 1 0 0;-25 0 0 0;0 0 0 1;0 0 -144 0];
gamma_e=.0001*ones(length(Am),1);
Fu=mats;
K_u=ones(size(Fu,1),1);
Theta=ones(1,size(Fu,1));
Bbar=padarray(B,[size(Fu,1),0],0,'post');

Abar=[Am B*Theta; zeros(size(Fu,1),size(Am,2)) Fu];
eigAbar=eig(Abar);
Cbar=padarray(C,[0,size(Fu,1)],0,'post');
rank(obsv(Abar,Cbar))
[Klqr,Slqr,elqr]=lqr(Abar',Cbar',1*eye(86),0.1,0)
Klqr=Klqr.'
[eig(Abar-Klqr*Cbar)]
K1=Klqr(1:length(Am),:);
K2=Klqr(length(Am)+1:end,:);
%% Plots and Sim
sim2=sim('real_UIO_EEG_2021_03_03_01.slx',180);
figure
ax1=subplot(1,3,1);
plot(sim2.y1)
hold on
plot(sim2.yhat1)
grid on
xlim([18,22])
legend('Actual','Predicted')
title('($y_1$ and $\hat{y}_1$)')
xlabel('Time (s)')
ylabel('Output State ($y$)')
ax2=subplot(1,3,2);
plot(sim2.y2)
hold on
plot(sim2.yhat2)
grid on
xlim([18,22])
legend('Actual','Predicted')
title('($y_2$ and $\hat{y}_2$)')
xlabel('Time (s)')
ylabel('Output State ($y$)')
ax3=subplot(1,3,3);
plot(sim2.ey.Time,sim2.ey.Data(:,5)) %SINGLE CHannel
grid on
xlim([18,22])
title('Error ($e_y$)')
xlabel('Time (s)')
ylabel('Observer Error ($e_y$)')
linkaxes([ax1,ax2,ax3],'y');
sgtitle('Comparison of Single Channels and Error')


figure
subplot(2,2,1)
plot(sim2.xhat1)
grid on
title('Estimated Internal State ($\hat{x}_1$)')
xlabel('Time (s)')
ylabel('Internal State ($\hat{x}_1$)')
subplot(2,2,2)
grid on
plot(sim2.xhat2)
title('Estimated Internal State ($\hat{x}_2$)')
xlabel('Time (s)')
ylabel('Internal State ($\hat{x}_2$)')
subplot(2,2,3)
plot(sim2.xhat3)
grid on
title('Estimated Internal State ($\hat{x}_3$)')
xlabel('Time (s)')
ylabel('Internal State ($\hat{x}_3$)')
subplot(2,2,4)
plot(sim2.xhat4)
grid on
title('Estimated Internal State ($\hat{x}_4$)')
xlabel('Time (s)')
ylabel('Internal State ($\hat{x}_4$)')

figure
subplot(1,1,1)
plot(sim2.uhat)
grid on
title('Estimated Unknown Input ($\hat{u}$)')
xlabel('Time (s)')
ylabel('Unknown Input ($\hat{u}$)')

