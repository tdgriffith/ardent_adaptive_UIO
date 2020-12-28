opengl software
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
A=[0 1; 0 0];
B=[0; 1];
C= [1 0];
Gamma=[1;1];
Theta=[1 0];
F=[0 1;-1 0];
Ahat=[A Gamma*Theta; zeros(size(A)), F];
Bhat=[B;0;0];
Chat=[C,0,0];
Khat=[4;3;2;-2];
G=[-1 -2];
Gd=[-3 -1];
Ghat=[G Gd];

out=sim('CLE3_sim.slx');

figure
subplot(2,2,1)
plot(out.x)
hold on
plot(out.xhat)
legend('$x_1$','$x_2$','$\hat{x}_1$','$\hat{x}_2$')
grid on
title('Comparison of Plant States ($x$ vs. $\hat{x}$)')
subplot(2,2,2)
plot(out.z)
hold 
plot(out.zhat)
legend('$z_1$','$z_2$','$\hat{z}_1$','$\hat{z}_2$')
grid on
title('Comparison of Distrubance States ($z$ vs. $\hat{z}$)')
subplot(2,2,3)
plot(out.y)
hold on
plot(out.yhat)
legend('$y$','$\hat{y}$')
grid on
title('Comparison of Output State ($y$ vs. $\hat{y}$)')
subplot(2,2,4)
plot(out.e)
legend('$e_{x1}$','$e_{x2}$','$e_{z1}$','$e_{z2}$')
grid on
title('Comparison of Error States')