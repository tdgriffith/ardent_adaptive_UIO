%% Unknown Input Observer with Adaptive State Estimation
% Implementation of our method:
%% Setup
%opengl software
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%% Prepare Data
fs=128;
dt=inv(fs);
EEG=load('F:\Documents\MATLAB\DEAP\s06.mat');
Y_data=squeeze(EEG.data([1],[2,8,14,20,25,30],:));

% t=EEG.samples([100:end],1);
% Y_data=detrend(Y_data);
% Y_data=bandpass(Y_data,[1 49],fs); % PLAY with this
% Y_data=Y_data';

    
Y_sim=array2timetable(Y_data(:,1:128*60)','SampleRate',fs);
Y_sim1=Y_sim(:,1);
Y_sim2=Y_sim(:,2);
Y_sim3=Y_sim(:,3);
Y_sim4=Y_sim(:,4);
Y_sim5=Y_sim(:,5);
Y_sim6=Y_sim(:,6);
%% OMA
order = 40;
s = 2*order;
opt_order=30;
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
%% 
Am=A_data{opt_order};
C=C_data{opt_order};
B=B2;


% Amc=[real(Am) -imag(Am);imag(Am) real(Am)];
% Bc=[real(B) -imag(B);imag(B) real(B)];
% Ac=[real(A) -imag(A);imag(A) real(A)];
gamma_e=[0.0001;0.0001;0.0001;0.0001;0.0001;0.0001]';
ku=-0.01; %Placement needs work
%ku=0
K_u=[-1;-10;-10;-100;-10];
K_u=[K_u,K_u,K_u,K_u,K_u,K_u];
K_u=-rand(5,6)

for i=1:1:100000
    K_u=-rand(5,6)*randi([-10,-1]);
    Ac_bar=[Am B*Theta; K_u*C Fu];
    eigs=eig(Ac_bar)
    if all(eig(Ac_bar)<0)
        K_ufinal=K_u
    else
    end
end
    

Abar=[Am B*Theta; zeros(size(Fu,1),size(Am,2)), Fu]
% Kbar=[zeros(size)]

Fu=[0 0 0 0 0; 0 0 1 0 0;0 -25 0 0 0;0 0 0 0 1;0 0 0 -144 0];
Theta=[1,1,1,1,1];
Ac_bar=[Am B*Theta; K_u*C Fu];
eig(Ac_bar)
%%
sim_out=sim('complex_UIO_alternative_2021_01_27.slx',60);