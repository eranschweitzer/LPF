clear variables; close all;
v = 0.7:0.01:1.3;
u = log(v);
y = 0:0.1:4;
[U,Y] = meshgrid(u,y);

Z = exp(U).*Y;
% Z = U.^2.*Y;
Z2 = (1 + U).*Y;
% Z2 = (2*U - 1).*Y;
figure;
mesh(U,Y,Z,zeros(size(Z)))
hold on;
surf(U,Y,Z2,ones(size(Z)))
xlabel('u')
ylabel('|y|')