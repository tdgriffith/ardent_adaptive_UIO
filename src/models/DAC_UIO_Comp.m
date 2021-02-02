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
A=A-0.1*eye(2);
B=[0;1];
C=[1, epsilon];
Phi_k = [1 cos(t) sin(t)].';
gamma_e=1;
gamma_u=2;
g1=10;
g2=2;
g3=-5;
ThetaLstar=[g1,g2,g3];
Amc=[real(Am) -imag(Am);imag(Am) real(Am)];
Bc=[real(B) -imag(B);imag(B) real(B)];
Ac=[real(A) -imag(A);imag(A) real(A)];
%% DAC with known A,B,C
sigma_e=0.1;
sigma_d=0.2;
%% Sim
sim_out=sim('UIOwDAC_balas_01.slx',60);
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

ax5=subplot(2,3,5);
plot(sim_out.ThetaL)
grid on

title('Basis $\alpha$')
xlabel('Time (s)')
ylabel('$\alpha$')
legend('$\alpha_{1}$','$\alpha_{\sin}$','$\alpha_{\cos}$')

%linkaxes([ax1,ax2],'y');
sgtitle(['UIO for Input of ', num2str(g2),'$\sin{t}+$',num2str(g3),'$\cos{t}+$',num2str(g1)])

%% DAC Plots
figure
ax6=subplot(2,2,1);
plot(sim_out.y_DAC)
grid on
%xlim([18,22])
legend('$y$')
title('Output State')
xlabel('Time (s)')
ylabel('Output State ($y$)')

ax7=subplot(2,2,2);
plot(sim_out.eu_DAC)
grid on
legend('$e_u$')
title('Rejected Input Error')
xlabel('Time (s)')
ylabel('Rejected Input Error ($e_u$)')

ax8=subplot(2,2,3);
plot(-1*sim_out.ThetaL_DAC)
grid on
legend('$\Theta_{const}$,$\Theta_{\sin}$,$\Theta_{cos}$')
title('Rejected Input Basis Amplitudes')
xlabel('Time (s)')
ylabel('Input Amplitudes ($e_u$)')

ax9=subplot(2,2,4);
plot(-1*sim_out.uhat_DAC)
hold on
plot(sim_out.u)
grid on
legend('$\hat{u}$')
title('Input Comparison')
xlabel('Time (s)')
ylabel('Input ($u$,$\hat{u}$)')

sgtitle(['Adaptive DAC for Input of ', num2str(g2),'$\sin{t}+$',num2str(g3),'$\cos{t}+$',num2str(g1)])
%% Quantum
P1=eye(2);
P2=[0 1;1 0];
P3=[1 0;0 -1];
P4=[0 -1j;1j 0];

A=0.2*P1+2*P2+1.5*P3+1*P4;
%A=0.2*P1;
A=A*-1j;
A=A-1*eye(2)
syms ep Lx
B=[0;1];
C_sym=[1, ep];
Am_sys=A-B*Lx*C_sym;
Am=double(subs(Am_sys,{ep,Lx},{0.1,1}))
%Am(1,1)=Am(1,1)-0.1;
C=double(subs(C_sym,{ep},{0.1}));
Amc=[real(Am) -imag(Am);imag(Am) real(Am)];
Bc=[real(B) -imag(B);imag(B) real(B)];
Ac=[real(A) -imag(A);imag(A) real(A)];

gamma_e=4;
gamma_u=.01;
g1=2;
g2=0+0j;
g3=-0+0j;
ThetaLstar=[g1,g2,g3];
%% Other algortihm
Am=[0 1;-1 -1];
C=[1 0.1];
B=[0;1];
Lx=0.2;
A=Am+B*Lx*C;

Amc=[real(Am) -imag(Am);imag(Am) real(Am)];
Bc=[real(B) -imag(B);imag(B) real(B)];
Ac=[real(A) -imag(A);imag(A) real(A)];
gamma_e=.0001;
K_u=-0.1;
g1=-1+3j;
g2=0;
g3=0;




