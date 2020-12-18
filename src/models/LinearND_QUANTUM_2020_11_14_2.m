% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

clear all, close all
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
% figpath = '../figures/';
% addpath('./utils');

%% generate Data
n=2;
basis=gen_iH_basis(n);
A = zeros(n,n);
for i=1:length(basis)
    A=A+(basis{i}*rand());
end
%A=[.25j 0.8+0.2j;-0.8+0.2j 0.75j];
%A=[-1, 0;0 0];
%A=[0, j;-j 0];
%A=A*j;

rhs = @(x)A*x;   % ODE right hand side
dt=.01;
tspan=[0:dt:10];   % time span
x0=rand(n,1)+j*rand(n,1); %initial conditions
%x0=[1j,-1j];
%x0=[-1;0]
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
[t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate


%% compute Derivative 
eps = 0.0;      % noise strength
dx=gradient(x.',dt).';
dx=dx+eps*randn(size(dx));

%% Custom Library
Theta_cust=[];
basis_mat=cell2mat(permute(basis,[1,3,2]));
basis_mat_sparse=[];
for i=1:length(basis_mat)
    basis_mat_sparse{i}=sparse(basis_mat(:,:,i));
end

Theta_cust=[];
% for i = 1:length(basis)
%     x_sigma = (basis{i})*B_ts.';
%     Theta_cust = [Theta_cust,x_sigma.'];
%     disp(i/length(basis))
% end
display('Theta Start');
tic
Theta_cust = reshape(reshape(permute(basis_mat,[2,1,3]),n,[]).'*x.',[],n);
%Theta_cust=Theta_cust.';
Theta_cust=reshape(Theta_cust,n^3,[]).';
Theta_cust=reshape(Theta_cust,[],length(basis));
Theta=[real(Theta_cust);imag(Theta_cust)];
toc
display('Theta End');

%% Estimate State Matrix
display('Least Sq Start');
tic
dx3=[real(dx(:));imag(dx(:))];
% G=gpuArray(Theta);
% b=gpuArray(dx3);
G=Theta;
b=dx3;
Xi=lsqr(G,b);
%Xi=gather(x);

lambda=0.0;
smallinds = (abs(Xi)<lambda);   % find small coefficients
biginds = ~smallinds;
G=G(:,biginds);
Xi=lsqr(G,b);
basis(smallinds) = [];


toc
display('Least Sq End');
%%

%n_offdiag=(((n*n)-n)/2)-1;
%Xi(end-n_offdiag:end)=Xi(end-n_offdiag:end)*-1;

A_solve=zeros(n,n);
for i = 1:length(basis)
    A_solve=A_solve+(Xi(i)*basis{i});
end
A_solve
A
error=norm(A-A_solve)

[tA,xA]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
rhs2 = @(x)A_solve*x;   % ODE right hand side
[tB,xB]=ode45(@(t,x)rhs2(x),tspan,x0,options);  % integrate
normx=norm(xA);
figure
subplot(1,2,1)
plot(tA,xA,'LineWidth',1.5)
hold on
plot(tB,xB,'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, $x_k$')
title('Real. Comp. A')

subplot(1,2,2)

plot(tA,imag(xA),'LineWidth',1.5)
hold on
plot(tB,imag(xB),'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, $x_k$')
title('Imag. Comp. A')
sgtitle(['Simulation Results True A vs. SysID A: Error ',num2str(error)])

%% Complex Sim
realA=real(A);
imagA=imag(A);
Ac=[realA -imagA;imagA realA];
realIC=real(x0);
imagIC=imag(x0);
simout=sim('complex_test.slx',10);

figure
subplot(1,2,1)
plot(tA,xA,'LineWidth',1.5)
hold on
plot(simout.realx1,'k--','LineWidth',1.2)
plot(simout.realx2,'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, $x_k$')
title('Real. Comp. A')

subplot(1,2,2)

plot(tA,imag(xA),'LineWidth',1.5)
hold on
plot(simout.imagx1,'k--','LineWidth',1.2)
plot(simout.imagx2,'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, $x_k$')
title('Imag. Comp. A')
sgtitle(['Simulation Results True A vs. SysID A: Error ',num2str(error)])

%% UIO Complex

