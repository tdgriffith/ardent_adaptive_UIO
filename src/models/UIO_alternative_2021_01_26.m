%% Setup
%opengl software
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%% Step Only
Am=[0 1;-1 -1];
C=[1 0.1];
B=[0;1];
Lx=0.2;
A=Am+B*Lx*C;

Amc=[real(Am) -imag(Am);imag(Am) real(Am)];
Bc=[real(B) -imag(B);imag(B) real(B)];
Ac=[real(A) -imag(A);imag(A) real(A)];
gamma_e=.0001;
K_u=-0.05;
g1=-1-2j;
g2=0; %STEP ONLY
g3=0; %STEP ONLY

%% Plots for step only
sim_out=sim('complex_UIO_3.slx',60);
figure
ax1=subplot(2,3,1);
plot(sim_out.y.Time,abs(squeeze(sim_out.y.Data)))
hold on
plot(sim_out.yhat.Time,abs(squeeze(sim_out.yhat.Data)))
grid on
%xlim([18,22])
legend('$y$','$\hat{y}$')
title('($y$ and $\hat{y}$)')
xlabel('Time (s)')
ylabel('Output State ($y$)')
ax2=subplot(2,3,2);
plot(sim_out.x.Time,abs(squeeze(sim_out.x.Data)))
hold on
plot(sim_out.xhat.Time,abs(squeeze(sim_out.xhat.Data)),'-.')
grid on
%xlim([18,22])
%legend('$\textit{R}(x_1)$','$\textit{R}(x_2)$','$\textit{R}(x_3)$','$\textit{R}(x_4)$','$\textit{R}( \hat{x}_1)$','$\textit{R}( \hat{x}_2)$','$\textit{R}( \hat{x}_3)$','$\textit{R}( \hat{x}_4)$','$\textit{C}(x_1)$','$\textit{C}(x_2)$','$\textit{C}(x_3)$','$\textit{C}(x_4)$','$\textit{C}( \hat{x}_1)$','$\textit{C}( \hat{x}_2)$','$\textit{C}( \hat{x}_3)$','$\textit{C}( \hat{x}_4)$')
title('($x$ and $\hat{x}$)')
xlabel('Time (s)')
ylabel('Internal State ($x$)')

ax3=subplot(2,3,3);
%plot(sim_out.ey)
plot(sim_out.ey.Time,abs(squeeze(sim_out.ey.Data)))
grid on
legend('$e_y$')
title('State Error')
xlabel('Time (s)')
ylabel('Error ($e_y$)')

ax4=subplot(2,3,4);
plot(sim_out.u.Time,abs(squeeze(sim_out.u.Data)))
hold on
plot(sim_out.u.Time,abs(squeeze(sim_out.uhat.Data)),'-.')
plot(sim_out.uhat_tot.Time,abs(squeeze(sim_out.uhat_tot.Data)),'-.')
grid on
title('Input State ($u$ and $\hat{u}$)')
xlabel('Time (s)')
ylabel('Input State')
%legend('$u$','$\hat{u}$')

ax5=subplot(2,3,5);
plot(sim_out.u.Time,abs(squeeze(sim_out.eu.Data)))
hold on
plot(sim_out.u.Time,abs(squeeze(sim_out.eu_tot.Data)))
grid on

title('Input Error ($e_u$)')
xlabel('Time (s)')
ylabel('Input Error ($e_u$)')
%legend('$e_u$')

ax6=subplot(2,3,6);
plot(sim_out.alpha)
grid on
hold on
plot(sim_out.alpha.Time,imag(squeeze(sim_out.alpha.Data)))
title('Basis $\alpha$')
xlabel('Time (s)')
ylabel('$\alpha$')
legend('$\alpha_{1_{re}}$','$\alpha_{1_{im}}$')

%linkaxes([ax1,ax2],'y');
sgtitle(['UIO for Input of ', num2str(g2),'$\sin{t}+$',num2str(g3),'$\cos{t}+$',num2str(g1)])

%% Cosine and Sine and Step
Am=[0 1;-1 -1];
C=[1 0.1];
B=[0;1];
Lx=0.2;
A=Am+B*Lx*C;
%A=Am

Amc=[real(Am) -imag(Am);imag(Am) real(Am)];
Bc=[real(B) -imag(B);imag(B) real(B)];
Ac=[real(A) -imag(A);imag(A) real(A)];
gamma_e=0.001;
ku=-0.2; %Placement needs work
K_u=[-.1;-.2;9];
g1=-5;
g2=2; 
g3=3; 
Fu=[0 0 0; 0 0 1 ;0 -25 0];
Theta=[1,1,1];
Ac_bar=[Am B*Theta; K_u*C Fu];
eig(Ac_bar)



%% 3 Signal Plots
sim_out=sim('complex_UIO_alternative_2021_01_27.slx',60);
figure
ax1=subplot(2,3,1);
plot(sim_out.y.Time,abs(squeeze(sim_out.y.Data)))
hold on
plot(sim_out.yhat.Time,abs(squeeze(sim_out.yhat.Data)))
grid on
%xlim([18,22])
legend('$y$','$\hat{y}$')
title('($y$ and $\hat{y}$)')
xlabel('Time (s)')
ylabel('Output State ($y$)')
ax2=subplot(2,3,2);
plot(sim_out.x.Time,abs(squeeze(sim_out.x.Data)))
hold on
plot(sim_out.xhat.Time,abs(squeeze(sim_out.xhat.Data)),'-.')
grid on
%xlim([18,22])
%legend('$\textit{R}(x_1)$','$\textit{R}(x_2)$','$\textit{R}(x_3)$','$\textit{R}(x_4)$','$\textit{R}( \hat{x}_1)$','$\textit{R}( \hat{x}_2)$','$\textit{R}( \hat{x}_3)$','$\textit{R}( \hat{x}_4)$','$\textit{C}(x_1)$','$\textit{C}(x_2)$','$\textit{C}(x_3)$','$\textit{C}(x_4)$','$\textit{C}( \hat{x}_1)$','$\textit{C}( \hat{x}_2)$','$\textit{C}( \hat{x}_3)$','$\textit{C}( \hat{x}_4)$')
title('($x$ and $\hat{x}$)')
xlabel('Time (s)')
ylabel('Internal State ($x$)')

ax3=subplot(2,3,3);
%plot(sim_out.ey)
plot(sim_out.ey.Time,abs(squeeze(sim_out.ey.Data)))
grid on
legend('$e_y$')
title('State Error')
xlabel('Time (s)')
ylabel('Error ($e_y$)')

ax4=subplot(2,3,4);
plot(sim_out.u.Time,abs(squeeze(sim_out.u.Data)))
hold on
plot(sim_out.u.Time,abs(squeeze(sim_out.uhat.Data)),'-.')
plot(sim_out.uhat_tot.Time,abs(squeeze(sim_out.uhat_tot.Data)),'-.')
grid on
title('Input State ($u$ and $\hat{u}$)')
xlabel('Time (s)')
ylabel('Input State')
%legend('$u$','$\hat{u}$')

ax5=subplot(2,3,5);
plot(sim_out.u.Time,abs(squeeze(sim_out.eu.Data)))
hold on
plot(sim_out.u.Time,abs(squeeze(sim_out.eu_tot.Data)))
grid on

title('Input Error ($e_u$)')
xlabel('Time (s)')
ylabel('Input Error ($e_u$)')
%legend('$e_u$')

%ax6=subplot(2,3,6);

sgtitle(['UIO Overview: Input of ', num2str(g2),'$\sin{t}+$',num2str(g3),'$\cos{t}+$',num2str(g1)])

figure
ax7=subplot(2,3,1)
plot(sim_out.u.Time,real(squeeze(sim_out.u.Data)))
hold on
plot(sim_out.uhat_tot.Time,real(squeeze(sim_out.uhat_tot.Data)))
grid on
legend('$u$','$\hat{u}+Ly$')
title('Total Input in Real ($\hat{u}+Ly$)')
xlabel('Time (s)')
ylabel('Input')

ax8=subplot(2,3,2)
plot(sim_out.u.Time,imag(squeeze(sim_out.u.Data)))
hold on
plot(sim_out.uhat_tot.Time,imag(squeeze(sim_out.uhat_tot.Data)))
grid on
legend('$u$','$\hat{u}+Ly$')
title('Total Input in Complex ($\hat{u}+Ly$))')
xlabel('Time (s)')
ylabel('Input')

ax9=subplot(2,3,3)
plot(sim_out.eu_tot.Time,abs(squeeze(sim_out.eu_tot.Data)))
grid on
legend('$u$','$\hat{u}+Ly$')
title('Abs Error ($u-(\hat{u}+Ly)$')
xlabel('Time (s)')
ylabel('Input Error')

ax10=subplot(2,3,4)
plot(sim_out.u.Time,real(squeeze(sim_out.u.Data)))
hold on
plot(sim_out.uhat_tot.Time,real(squeeze(sim_out.uhat.Data)))
grid on
legend('$u$','$\hat{u}$')
title('Estimated Input in Real ($\hat{u}$)')
xlabel('Time (s)')
ylabel('Input')

ax11=subplot(2,3,5)
plot(sim_out.u.Time,imag(squeeze(sim_out.u.Data)))
hold on
plot(sim_out.uhat_tot.Time,imag(squeeze(sim_out.uhat.Data)))
grid on
legend('$u$','$\hat{u}$')
title('Estimated Input in Complex ($\hat{u}$)')
xlabel('Time (s)')
ylabel('Input')

ax12=subplot(2,3,6)
plot(sim_out.eu_tot.Time,abs(squeeze(sim_out.eu.Data)))
grid on
legend('$u$','$\hat{u}$')
title('Abs Error ($u-(\hat{u}))$')
xlabel('Time (s)')
ylabel('Input Error')
sgtitle(['UIO Input Details: ', num2str(g2),' $\sin{t}+$',num2str(g3),' $\cos{t}+$',num2str(g1)])
%% 5 unknown inputs
Am=[0 1;-1 -1];
C=[1 0.1];
B=[0;1];
Lx=0.2;
A=Am+B*Lx*C;
A=Am

Amc=[real(Am) -imag(Am);imag(Am) real(Am)];
Bc=[real(B) -imag(B);imag(B) real(B)];
Ac=[real(A) -imag(A);imag(A) real(A)];
gamma_e=0.0001;
ku=-0.1; %Placement needs work
%ku=0
K_u=[ku;ku-.01;ku;ku-.1;ku-.1];
g1=-1+2j;
g2=2+3j; 
g3=4+1j; 
g4=20
g5=-10j
Fu=[0 0 0 0 0; 0 0 1 0 0;0 -25 0 0 0;0 0 0 0 1;0 0 0 -144 0];
Theta=[1,1,1,1,1];
Ac_bar=[Am B*Theta; K_u*C Fu];
eig(Ac_bar)


%%
% syms ku1 ku2 ku3
% Ku=[ku1;ku2;ku3];
% a0=-1;
% a1=-1;
% a2=-26+ku3/10+ku2/10+ku1/10;
% a3=-3/2*ku2+11/10*ku3+ku1-25;
% a4=5/2*ku1-25*ku2+ku3-25;
% a5=25*ku1;
% b1=(a1*a2-a0*a3)/a1;
% b2=(a1*a4-a0*a5)/a1;
% b3=(a1*0-a0*0)/a1;
% c1=(b1*a3-a1*b2)/b1;
% c2=(b1*a5-a1*b3)/b1;
% c3=(b1*0-a1*0)/b1;
% d1=(c1*b2-b1*c2)/c1;
% d2=(c1*b3-b1*c3)/c1;
% d3=(c1*0-0*0)/c1;
% e1=(d1*c2-c1*d2)/d1;
% e2=(d1*c3-c1*d3)/d1;
% e3=(d1*0-c1*0)/d1;
% routh=[a0 a2 a4;a1 a3 a5;b1 b2 b3;c1 c2 c3;d1 d2 d3;e1 e2 e3]
%%
% Abar=[Am,B*Theta;zeros(3,2), Fu];
% Cbar=[C,0,0,0];
% cvx_begin sdp
%     variable P(5,5) symmetric
%     variable K(3,1)
%     Abar'*P+P*Abar+Cbar'*(P*[0;0;K])'+(P*[0;0;K])*Cbar <= -1e-5*eye(5)
%     P >= 1e-5*eye(5)
% cvx_end
% 
% K=inv(P)*[0;0;Y];
% mat=Abar'*P+P*Abar+Cbar'*Y'+Y*Cbar;
% try chol(-mat)
%     disp('Matrix is symmetric positive definite.')
% catch ME
%     disp('Matrix is not symmetric positive definite')
% end
%% Parametric poles ew
ku=-9/25; %Placement needs work
%ku=0
K_u=[ku;ku;ku];
Fu=[0 0 0; 0 0 1 ;0 -25 0];
Theta=[1,1,1];
Ac_bar=[Am B*Theta; K_u*C Fu];
eig(Ac_bar)
v_ts=[];
ku_ts=[];
for ku=-15:0.1:15
    ku_ts=[ku_ts,ku]
    K_u=[-.1;-.2;9]
    Ac_bar=[Am B*Theta; K_u*C Fu]
    v=eig(Ac_bar)
    v_ts=[v_ts,v]
end
figure
plot(ku_ts,real(v_ts))
