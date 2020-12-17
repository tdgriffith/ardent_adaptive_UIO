%% Unknown Input Observer with Adaptive State Estimation
% Implementation from the following paper:
% B. Alenezi, J. Hu and S. H. {ak, "Adaptive Unknown Input and State Observers," 2019 American Control Conference (ACC), Philadelphia, PA, USA, 2019, pp. 2434-2439, doi: 10.23919/ACC.2019.8815288.
%% Setup
%opengl software
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
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
%% DAC with known A,B,C

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


