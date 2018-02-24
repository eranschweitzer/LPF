function vt_comp_plots(vars, vtrue)
%%% assumption is that theta is given in radians!!!!
figure;
plot(vtrue.v,'-o','linewidth',3);
hold on;
plot(vars.v,'-*','linewidth',3)
legend('true', 'linpf')
xlabel('bus number')
ylabel('p.u. voltage')
ax = gca;
ax.FontSize = 16;

figure;
plot(vtrue.t*180/pi,'-o','linewidth',3);
hold on;
plot(vars.theta*180/pi,'-*','linewidth',3)
legend('true', 'linpf')
xlabel('bus number')
ylabel('bus angle [degree]')
ax = gca;
ax.FontSize = 16;