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
EEG=load('F:\Documents\MATLAB\group_01_mat\group_01_sujet_01.mat');
Y_data=EEG.samples([100:end],[2,8,16,20,24,30]);

t=EEG.samples([100:end],1);
Y_data=detrend(Y_data);
Y_data=bandpass(Y_data,[4 45],fs);
Y_data=Y_data';

    
Y_sim=array2timetable(Y_data(:,1:512*60)','SampleRate',fs);
Y_sim1=Y_sim(:,1);
Y_sim2=Y_sim(:,2);
Y_sim3=Y_sim(:,3);
Y_sim4=Y_sim(:,4);
Y_sim5=Y_sim(:,5);
Y_sim6=Y_sim(:,6);
%% OMA
order = 20;
s = 2*order;
opt_order=20;
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


%% Solve First LMIs
cvx_begin sdp %CONTINUOUS
    variable P(opt_order,opt_order) symmetric
    variable Y(size(C_data{opt_order},2),size(C_data{opt_order},1))
    A_data{opt_order}'*P+P*A_data{opt_order}-C_data{opt_order}'*Y'-Y*C_data{opt_order} <= -0.1*eye(opt_order);
    P >= 1*eye(opt_order);
cvx_end

%% DISCRETE LMI
cvx_begin sdp
    variable P(opt_order,opt_order) symmetric
    variable Y(size(C_data{opt_order},2),size(C_data{opt_order},1))
    [P, (P*A_data{opt_order}-Y*C_data{opt_order})' ;
        P*A_data{opt_order}-Y*C_data{opt_order}, P] >= 1e-16*eye(2*opt_order)
cvx_end
%%

L=inv(P)*Y;
%mat=A_data{opt_order}'*P+P*A_data{opt_order}-C_data{opt_order}'*Y'-Y*C_data{opt_order};
abs(eig(A_data{opt_order}-L*C_data{opt_order}))

F=B2'*P/(C_data{opt_order});

%%
Gamma=.001;
Abar=[A_data{opt_order}-L*C_data{opt_order}, B2;
    -Gamma*B2'*P zeros(1)];
abs(eig(Abar))



