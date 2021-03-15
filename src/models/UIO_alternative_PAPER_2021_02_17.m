%% Setup: Script to generate paper figures
%opengl software
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%% Cosine and Sine and Step
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
ku=-0.2; %Placement needs work
K_u=[-.1;-.2;9];
g1=1;
g2=2; 
g3=-3; 
Fu=[0 0 0; 0 0 1 ;0 -25 0];
Theta=[1,1,1];
Ac_bar=[Am B*Theta; K_u*C Fu];
eig(Ac_bar)
Abar=[Am B*Theta; zeros(size(Fu,1),size(Am,2)) Fu];
Cbar=[C,0,0,0];
rank(obsv(Abar,Cbar))
Bbar=[B;0;0;0];

% cvx_begin sdp
%     variable P(5,5) symmetric
%     variable Y(5,1)
%     minimize(trace(P));
%     subject to
%         Abar'*P+P*Abar-Cbar'*Y'-Y*Cbar <= -0.1*eye(5)
%         P >= 1*eye(5)
%         sum(P(1,:)) >= 1000
%         %abs(Y(1,1))<=0.1
%         %abs(Y(2,1))<=0.1
%         sum(P(2,:)) >= 1000
%         %P(2,2) >= 10
%         %P*Bbar==Cbar'
% cvx_end
cvx_begin sdp
    variable P(5,5) symmetric
    variable Y(5,1)
    minimize(norm(P));
    subject to
        Abar'*P+P*Abar-Cbar'*Y'-Y*Cbar <= -.001*eye(5)
        P >= 1*eye(5)
        sum(P(1,:)) >= 1000
        abs(Y(1,1))<=0.1
        abs(Y(2,1))<=0.1
        %sum(P(2,:)) >= 1000
        %P(2,2) >= 10
        %P*Bbar==Cbar'
cvx_end

% cvx_begin sdp
%     variable P(5,5) symmetric
%     variable Y(5,1)
%         Abar'*P+P*Abar-Cbar'*Y'-Y*Cbar <= -.001*eye(5)
%         P >= 1*eye(5)
%         sum(P(1,:)) >= 0
%         abs(Y(1,1))<=0.1
%         abs(Y(2,1))<=0.1
%         sum(P(2,:)) >= 1000
%         %P(2,2) >= 10
%         %P*Bbar==Cbar'
% cvx_end
Kcvx=inv(P)*Y;
Kcvxbar=[Kcvx];
Kcvx(1,1)=0;
Kcvx(2,1)=0;
Kcvxbar2=[Kcvx];
[eig(Abar-Kcvxbar*Cbar),eig(Abar-Kcvxbar2*Cbar)]
[Kcvxbar,Kcvxbar2]


%% 3 Signal Plots
sim_out=sim('real_UIO_alternative_2021_02_17.slx',60);
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
%% Cosine and Sine and Step
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
%K_u=[-1;-1;3;-10;-1.6;10];

g1=2;
g2=1; 
g3=-3; 
g4=-1;
Fu=[0 0 0 0 0 0; 0,0,0,0,0, 0;0 0 0 1 0 0;0 0 -25 0 0 0;0 0 0 0 0 1;0 0 0 0 -144 0];
K_u=ones(size(Fu,1),1);
Theta=[1,1,1,1,1,1];
Ac_bar=[Am B*Theta; K_u*C Fu];
Kbar=[0;0;K_u];
Bbar=[B;0;0;0;0];
eig(Ac_bar)
Abar=[Am B*Theta; zeros(size(Fu,1),size(Am,2)) Fu];
Cbar=[C,0,0,0,0,0,0];
rank(obsv(Abar,Cbar))

%% Pole Placement
[Actrb,Bctrb,Cctrb,Tctrb,kctrb] = ctrbf(Abar,Kbar*Cbar,Cbar);
Ac_ctrb=Actrb(4:8,4:8);
Bc_ctrb=Bctrb(4:8,:);
rank(ctrb(Ac_ctrb,Bc_ctrb))
eig(Ac_ctrb)
Kctrb=place(Ac_ctrb,Bc_ctrb,[-2+3j,-2-3j,-10+1j,-10-1j,-6]);
eig(Ac_ctrb-Bc_ctrb*Kctrb)
K_out=padarray(Kctrb,[0,3],0,'pre');
K_out=K_out*Tctrb
eig(Abar-Kbar*Cbar*K_out) %works but wrong shape

%% Pole placement 2
[Actrb,Bctrb,Cctrb,Tctrb,kctrb] = ctrbf(Abar',Cbar',eye(8));
Ac_ctrb=Actrb(2:8,2:8);
Bc_ctrb=Bctrb(2:8,:);
rank(ctrb(Ac_ctrb,Bc_ctrb))
eig(Ac_ctrb)
Kctrb=place(Ac_ctrb,Bc_ctrb,[-2+3j,-2-3j,-10+1j,-10-1j,-6,-20,-5]);
[Klqr,Slqr,elqr]=lqr(Ac_ctrb,Bc_ctrb,eye(7),1,0)
eig(Ac_ctrb-Bctrb(2:8,:)*Kctrb)
eig(Ac_ctrb-Bctrb(2:8,:)*Klqr)
K_out=padarray(Kctrb,[0,1],0,'pre');
K_out_lqr=padarray(Klqr,[0,1],0,'pre');
eig(Actrb-Bctrb*K_out)
eig(Actrb-Bctrb*K_out_lqr)
K_out2=K_out*Tctrb;
K_out2_lqr=(K_out_lqr*Tctrb)'
K_out2=K_out2'
[eig(Abar),eig(Abar-K_out2*Cbar),eig(Abar-K_out2_lqr*Cbar)] 
[K_out2,Kbar]

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
%K_u=[-1;-1;3;-10;-1.6;10];

g1=2;
g2=1; 
g3=-3; 
g4=-1;
Fu=[0,0,0,0, 0;0 0 1 0 0;0 -25 0 0 0;0 0 0 0 1;0 0 0 -144 0];
K_u=ones(size(Fu,1),1);
Theta=[1,1,1,1,1];
Ac_bar=[Am B*Theta; K_u*C Fu];
Kbar=[0;0;K_u];
Bbar=[B;0;0;0;0;0];
eig(Ac_bar)
Abar=[Am B*Theta; zeros(size(Fu,1),size(Am,2)) Fu];
Cbar=[C,0,0,0,0,0];
rank(obsv(Abar,Cbar))

[Actrb,Bctrb,Cctrb,Tctrb,kctrb] = ctrbf(Abar',Cbar',eye(7));
Ac_ctrb=Actrb(3:7,3:7);
Bc_ctrb=Bctrb(3:7,:);
rank(ctrb(Ac_ctrb,Bc_ctrb))
eig(Ac_ctrb)
Kctrb=place(Ac_ctrb,Bc_ctrb,[-2+3j,-2-3j]);
cost_mat=eye(5);
[Klqr,Slqr,elqr]=lqr(Ac_ctrb,Bc_ctrb,cost_mat,1,0)
eig(Ac_ctrb-Bc_ctrb*Kctrb)
eig(Ac_ctrb-Bc_ctrb*Klqr)
K_out=padarray(Kctrb,[0,5],0,'pre');
K_out_lqr=padarray(Klqr,[0,2],0,'pre');
eig(Actrb-Bctrb*K_out)
eig(Actrb-Bctrb*K_out_lqr)
K_out2=(K_out*Tctrb)';
K_out2_lqr=(K_out_lqr*Tctrb)'

[Klqr2,Slqr,elqr]=lqr(Abar',Cbar',cost_mat,1,0)
Kout3=Klqr2';
Kout3(1,1)=0;
Kout3(2,1)=0;
%[eig(Abar),eig(Abar-K_out2*Cbar),eig(Abar-K_out2_lqr*Cbar)] 
[eig(Abar-Kout3*Cbar),eig(Abar-Klqr2'*Cbar)]


%% TABLE 5.3-2 Optimal Feedback Gain pg 404 Aircraft cOntrol and sim
A=[-0.3220 0.0640 0.0364 -0.9917 0.0003 0.0008 0;
    0 0 1 0.0037 0 0 0;
    -30.6492 0 -3.6784 0.6646 -0.7333 0.1315 0;
    8.5396 0 -0.0254 -0.4764 -0.0319 -0.0620 0;
    0 0 0 0 -20.2 0 0;
    0 0 0 0 0 -20.2 0;
    0 0 0 57.2958 0 0 -1];
B=[0 0;0 0;0 0;0 0;20.2 0;0 20.2;0 0];
C=[0 0 0 57.2958 0 0 -1;
    0 0 57.2958 0 0 0 0;
    57.2958 0 0 0 0 0 0;
    0 57.2958 0 0 0 0 0];
Q=[50 0 0 0 0 0 0;
    0 100 0 0 0 0 0;
    0 0 100 0 0 0 0;
    0 0 0 50 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 1];
rho=0.1;
Kbook=[-0.56 -0.44 0.11 -0.35;
    -1.19 -0.21 -0.44 0.26];
K=[-1 -1 -1 -1;-1 -1 -1 -1];
%K=[0 -0.55 0 -0.49;-1.14 0 0.55 0]+0.3*rand
eig(A-B*K*C)
X=eye(7);
R=rho*eye(2);
Ac=A-B*K*C;

%%
e=1;
k=0;
J=0;
evec=[];
while e > 0.1
    cvx_begin sdp
        variable P(7,7) symmetric
        variable S(7,7) 
        Ac*P+P*Ac+C'*K'*R*K*C <= -Q;
        Ac*S+S*Ac <= -X;
    cvx_end
    deltaK=inv(R)*B'*P*S*C'*inv(C*S*C')-K;
    K=K+0.00001*deltaK;
    Ac=A-B*K*C;
    if ~all(eig(Ac)<0)
        break
        display('unstable')
    end
    e=abs(1/2*trace(P*X)-J)
    evec=[e,evec];
    J=1/2*trace(P*X);
    k=k+1
end
%%
e=1;
k=0;
J=0;
evec=[];
while e > 0.02
    cvx_begin sdp
        variable P(7,7) symmetric
        variable S(7,7) 
        Ac*P+P*Ac+C'*K'*R*K*C <= -Q;
        Ac*S+S*Ac <= -X;
    cvx_end
    deltaK=inv(R)*B'*P*S*C'*inv(C*S*C')-K;
    K=K+0.00001*deltaK;
    Ac=A-B*K*C;
    if ~all(eig(Ac)<0)
        break
        display('unstable')
    end
    e=abs((1/2*trace(P*X)+1000*K(1,1)+K(1,2)+1000*K(1,3)+K(1,4)+K(2,1)+1000*K(2,2)+K(2,3)+1000*K(2,4)))
    evec=[e,evec];
    J=(1/2*trace(P*X)+1000*K(1,1)+K(1,2)+1000*K(1,3)+K(1,4)+K(2,1)+1000*K(2,2)+K(2,3)+1000*K(2,4));
    k=k+1
end

%%
cvx_begin sdp
    variable P(7,7) symmetric
    variable Y(7,4)
    A'*P+P*A-C'*Y'-Y*C <= -Q
    P >= 1e-5*eye(7)
cvx_end

K=pinv(B)*inv(P)*Y;
eig(A-B*K*C)
%%
%%
cvx_begin sdp
    variable P(7,7) symmetric
    variable S(7,7) symmetric
    variable Y(7,4)
    variable Z(2,7)
    minimize(0.5*trace(P*X));
    subject to
        A'*P+P*A-C'*Y'-Y*C  == -2*eye(7)
        A*S-B*Z+S*A'-Z'*B' == -X
        P >= 0.1*eye(7)
        S >= 0.1*eye(7)
cvx_end

K=pinv(B)*inv(P)*Y;
eig(A-B*K*C)

%K_cvx=pinv(B)*inv(P)*Y
%eig(A-B*K_cvx*C)
%% Parametric poles ew
ku=-9/25; %Placement needs work
%ku=0
K_u=[ku;ku;ku;ku];
v_ts=[];
ku_ts=[];
for ku=-15:0.1:15
    ku_ts=[ku_ts,ku]
    K_u=[2;10;-1.6;10]
    Ac_bar=[Am B*Theta; K_u*C Fu]
    v=eig(Ac_bar)
    v_ts=[v_ts,v]
end
figure
plot(ku_ts,real(v_ts))

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

