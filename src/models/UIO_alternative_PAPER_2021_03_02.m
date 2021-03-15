%% Setup: Script to generate paper figures
%opengl software
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%% So this version works with error feedback to both plant and input estimator
Am=[0 1;-1 -1];
C=[1 0.1];
B=[0;1];
Lx=-0.5;
A=Am+B*Lx*C;
%A=Am

% Amc=[real(Am) -imag(Am);imag(Am) real(Am)];
% Bc=[real(B) -imag(B);imag(B) real(B)];
% Ac=[real(A) -imag(A);imag(A) real(A)];
gamma_e=0.0001;
%ku=-0.1; %Placement needs work
%K_I=[2;10;-1.6;10];

g1=2;
g2=4; 
g3=1; 
g4=-2;
g5=0;
Fu=[0 1 0 0;-25 0 0 0;0 0 0 1;0 0 -144 0];
K_u=ones(size(Fu,1),1);
Theta=[1,1,1,1];
Ac_bar=[Am B*Theta; K_u*C Fu];
%Kbar_I=[0;0;K_I];
Bbar=[B;0;0;0;0];
eig(Ac_bar)
Abar=[Am B*Theta; zeros(size(Fu,1),size(Am,2)) Fu];
Cbar=[C,0,0,0,0];
rank(obsv(Abar,Cbar))


%%
cost_mat=diag([0.1 1 0.8 0.01 4 8]);
[Klqr2,Slqr,elqr]=lqr(Abar',Cbar',cost_mat,1,0)
Kout3=Klqr2';
[eig(Abar-Kout3*Cbar),eig(Abar-Klqr2'*Cbar)]
K1=Kout3(1:length(Am),:);
K2=Kout3(length(Am)+1:end,:);

%% Paper plots: x1 x2 ey u eu
sim_out=sim('real_UIO_alternative_2021_02_17_04.slx',60);

figure
ax1=subplot(2,1,1);
plot(sim_out.x.Time,squeeze(sim_out.ex.Data(1,:)),'Color', 'k', 'LineStyle', '-','LineWidth',1.5)
grid on
xlabel('Time (s)')
ylabel('$e_{x,1}$')
legend('$e_{x,1}$')
ax2=subplot(2,1,2);
plot(sim_out.x.Time,squeeze(sim_out.ex.Data(2,:)),'Color', 'k', 'LineStyle', '-','LineWidth',1.5)
grid on
title('$e_{x,2}$')
xlabel('Time (s)')
ylabel('$e_{x,2}$')
legend('$e_{x,2}$')

figure
ax3=subplot(2,1,1);
plot(sim_out.u.Time,squeeze(sim_out.u.Data),'Color', 'r', 'LineStyle', '-','LineWidth',1.5)
grid on
hold on
plot(sim_out.u.Time,squeeze(sim_out.uhat.Data),'Color', 'k', 'LineStyle', '--','LineWidth',0.8)
%title('True and Estimated Input vs. Time')
xlabel('Time (s)')
ylabel('Input')
legend('$u$','$\hat{u}$')
ax4=subplot(2,1,2);
plot(sim_out.eu.Time,squeeze(sim_out.eu.Data),'Color', 'k', 'LineStyle', '-','LineWidth',1.5)
grid on
xlabel('Time (s)')
ylabel('Input Error')
legend('$e_u$')



