% r/x ratio investigation 
% the basic relation:
% P = b( (rx)/2*delta^2 + delta)

%% case 1 b>0
% note that this is the normal case

rx    = linspace(-2,2);
delta = linspace(-pi/2,pi/2);

[RX,DELTA] = meshgrid(rx,delta);

P = 0.5*RX.*DELTA.^2 + DELTA;

figure;
pcolor(RX,DELTA,double(P>0));
xlabel('r/x')
ylabel('\theta_f - \theta_t')

figure;
mesh(RX,DELTA,P);
xlabel('r/x')
ylabel('\theta_f - \theta_t')
hold on;
surf(RX,DELTA,zeros(size(RX)));
mesh(RX,DELTA,DELTA);

figure;
mesh(RX,DELTA, DELTA-P)
xlabel('r/x')
ylabel('\theta_f - \theta_t')
%%
rx    = linspace(0,1);
delta = linspace(-pi/2,pi/2);

[RX,DELTA] = meshgrid(rx,delta);

P = 0.5*RX.*DELTA.^2 + DELTA;

figure;
mesh(RX,DELTA,P);
xlabel('r/x')
ylabel('\theta_f - \theta_t')