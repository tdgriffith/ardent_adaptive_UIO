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
setlmis([])
P=lmivar(1,[4,1]);
Y = lmivar(2,[4 2]); 
lmiterm([1 1 1 P],A,1,'s')
lmiterm([1 1 1 Y],-1,C,'s')
%lmiterm([1 1 1 -Y],C',-1)
lmiterm([-2,1,1,P],1,1)
lmis=getlmis;
options = [0,0,40,0,0];
[tmin,xfeas] = feasp(lmis);
P=dec2mat(lmis,xfeas,P)
Y=dec2mat(lmis,xfeas,Y);
L=inv(P)*Y
try chol(-1*(A'*P+P*A-C'*Y'-Y*C))
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end
%% Smaller system
C=1;
L=1;
R1=1;
R2=1;
R0=0.1;
A=[(R0+R1)/-L, -1/L;1/C,-1/(C*R2)];
C=[R0, 1];
setlmis([])
P=lmivar(1,[2,1]);
Y = lmivar(2,[2 1]); 
lmiterm([1 1 1 P],A,1,'s')
lmiterm([1 1 1 Y],-1,C,'s')
lmiterm([1 1 1 -Y],C',-1)
lmiterm([-2,1,1,P],1,1)
lmis=getlmis;
options = [0,0,10,0,0];
[tmin,xfeas] = feasp(lmis,options,-1);
%[tmin,xfeas] = feasp(lmis)
P=dec2mat(lmis,xfeas,P);
Y=dec2mat(lmis,xfeas,Y);
L=inv(P)*Y
%% Bigger
A=[-2379.2 0 0 0 0;0 -2.3 0 0.21 0;
    0 0 -2.3 0 0.21; 0 267.5 0 -43.83 0;
    0 0 267.54 0 -43.83];
B2=[68245 0 0;0 -2 0;0 0 2;0 232.75 0; 0 0 -232.75];
C=[1 0 0 0 0;0 0 0 1 0; 0 0 0 0 1];
setlmis([])
P=lmivar(1,[5,1]);
Y = lmivar(2,[5 3]); 
lmiterm([1 1 1 P],A,1,'s')
lmiterm([1 1 1 Y],-1,C,'s')
%lmiterm([1 1 1 -Y],C',-1)
lmiterm([-2,1,1,P],1,1)
lmis=getlmis;
options = [0,0,10,0,0];
[tmin,xfeas] = feasp(lmis,options,-1);
%[tmin,xfeas] = feasp(lmis)
P=dec2mat(lmis,xfeas,P);
Y=dec2mat(lmis,xfeas,Y);
L=inv(P)*Y