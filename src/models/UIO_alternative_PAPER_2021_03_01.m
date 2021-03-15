%% Setup: Script to generate paper figures
%opengl software
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%% Pole placement 3: Two is working, but kout has nonzero top values...
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
ku=-0.1; %Placement needs work
K_I=[2;10;-1.6;10];

g1=2;
g2=0; 
g3=0; 
g4=0;
Fu=[0 1 0 0;-25 0 0 0;0 0 0 1;0 0 -144 0];
K_u=ones(size(Fu,1),1);
Theta=[1,1,1,1];
Ac_bar=[Am B*Theta; K_u*C Fu];
Kbar_I=[0;0;K_I];
Bbar=[B;0;0;0;0];
eig(Ac_bar)
Abar=[Am B*Theta; zeros(size(Fu,1),size(Am,2)) Fu];
Cbar=[C,0,0,0,0];
rank(obsv(Abar,Cbar))

% [Actrb,Bctrb,Cctrb,Tctrb,kctrb] = ctrbf(Abar',Cbar',eye(7));
% Ac_ctrb=Actrb(1:7,1:7);
% Bc_ctrb=Bctrb(1:7,:);
% rank(ctrb(Ac_ctrb,Bc_ctrb))
% eig(Ac_ctrb)
% Kctrb=place(Ac_ctrb,Bc_ctrb,[-2+3j,-2-3j]);
% cost_mat=eye(5);
% [Klqr,Slqr,elqr]=lqr(Ac_ctrb,Bc_ctrb,cost_mat,1,0)
% eig(Ac_ctrb-Bc_ctrb*Kctrb)
% eig(Ac_ctrb-Bc_ctrb*Klqr)
% K_out=padarray(Kctrb,[0,5],0,'pre');
% K_out_lqr=padarray(Klqr,[0,2],0,'pre');
% eig(Actrb-Bctrb*K_out)
% eig(Actrb-Bctrb*K_out_lqr)
% K_out2=(K_out*Tctrb)';
% K_out2_lqr=(K_out_lqr*Tctrb)'
%%
cost_mat=diag([0.1 0.1 .01 .001 0.1 .001]);
[Klqr2,Slqr,elqr]=lqr(Abar',Cbar',cost_mat,1,0)
Kout3=Klqr2';
Kout3(1,1)=0;
Kout3(2,1)=0;
%[eig(Abar),eig(Abar-K_out2*Cbar),eig(Abar-K_out2_lqr*Cbar)] 
[eig(Abar+Kbar_I*Cbar),eig(Abar-Kout3*Cbar),eig(Abar-Klqr2'*Cbar)]
%%
cost_mat=diag([0.1 0.1 .01 .001 0.1 .001])
cvx_begin sdp
    variable P(6,6) symmetric
    variable Y(6,1)
    minimize(norm(P(1:6,1:6)));
    subject to
        Abar'*P+P*Abar-Cbar'*Y'-Y*Cbar <= -0.01*eye(6)
        P >= 0.01*eye(6)
        %P(1,1) >= 100
        %P(2,2) >= 10
        abs(Y(1,1))<=0.01
        abs(Y(2,1))<=0.01
        sum(P(2,:)) >= 1000
        %P(2,2) >= 10
        %P*Bbar==Cbar'
cvx_end
%%
cvx_begin sdp
    variable P(6,6) symmetric
    variable Y(6,1)
    Abar'*P+P*Abar-Cbar'*Y'-Y*Cbar <= diag([-.01 -.01 -.01 -1 -1 -1])
    P >= 0.1*eye(6)
    P(1,1) >= 0.1
    P(2,2) >= 100
    abs(Y(1,1))<=.1
    abs(Y(2,1))<=.1
    %sum(P(2,:)) >= 10
    %P*Bbar==Cbar'
cvx_end
K_cvx=inv(P)*Y;
K_cvx2=[0;0;K_cvx(3:end,:)];
[[0;0;K_I],K_cvx,K_cvx2]
[eig(Abar+Kbar_I*Cbar),eig(Abar-K_cvx*Cbar),eig(Abar-K_cvx2*Cbar)]
K_u=-[K_cvx(3:end,:)];

%% 3 Signal Plots
sim_out=sim('real_UIO_alternative_2021_02_17_02.slx',60);
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

%% Paper plots: x1 x2 ey u eu
sim_out=sim('real_UIO_alternative_2021_02_17_02.slx',60);

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

%% Toy example
A=[4 0 0;-4 0 -2; -2 -2 0];
B=[1;0;-1];
C=[1 4 16;0 -2 -12;-1 -2 -4];
[Abar,Bbar,Cbar,T,k] = ctrbf(A,B,C);
Ac=Abar(2:3,2:3);
Bc=Bbar(2:3,:);
Kctrb=place(Ac,Bc,[-1,-3]);
eig(Ac-Bc*Kctrb)
Kctrb_out=[0,Kctrb]*T;
eig(A-B*Kctrb_out)

