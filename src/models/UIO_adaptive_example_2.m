%% Unknown Input Observer with Adaptive State Estimation
% Implementation from the following paper:
% B. Alenezi, J. Hu and S. H. {ak, "Adaptive Unknown Input and State Observers," 2019 American Control Conference (ACC), Philadelphia, PA, USA, 2019, pp. 2434-2439, doi: 10.23919/ACC.2019.8815288.
%% Setup
%opengl software
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%% Define system
Ta=0.36;
alpha=1;
sigma2=0.00029;
wb2=11.77^2;
zetaa=0.6;
kb=-9.91;
wa=180;
kq=-0.07;

A=[-1/Ta, (alpha+sigma2*wb2)/Ta, (-kb*sigma2*wb2)/Ta, -kb*sigma2*wb2;
    -(1+wb2*Ta^2)/(Ta*(1+sigma2*wb2)), 1/Ta, ((Ta^2-sigma2)*kb*wb2)/(Ta*(1+sigma2*wb2)), 0;
    0, 0, 0, 1;
    0, 0, -wa^2, -2*zetaa*wa];
B1=[0;0;0;kq*wa^2];
B2=[1;0;1;0];
C=[1 0 0 0;0 1 0 0];
%% Solve First LMIs
cvx_begin sdp
    variable P(4,4) symmetric
    variable Y(4,2)
    A'*P+P*A-C'*Y'-Y*C <= -1e-5*eye(4)
    P >= 1e-5*eye(4)
cvx_end

L=inv(P)*Y;
mat=A'*P+P*A-C'*Y'-Y*C;
try chol(-mat)
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end

F=B2'*P/(C);

%%
Gamma=20;
Abar=[A-L*C, B2;
    -Gamma*B2'*P zeros(1)];

cvx_begin sdp
    variable P(4,4) symmetric
    variable Y(4,2)
    minimize(trace(P)+abs(100*Y(2,2)));
    subject to
        A'*P+P*A-C'*Y'-Y*C <= -1e-5*eye(4)
        P >= 1e-5*eye(4)
cvx_end