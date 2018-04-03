function vt_comp_plots(vars, vtrue, E, residual,bidx)
%%% assumption is that theta is given in radians!!!!
figure;
subplot(4,1,1)
plot(vtrue.v,'-o','linewidth',3);
hold on;
plot(vars.v,'-*','linewidth',3)
legend('true', 'linpf')
xlabel('bus number')
ylabel('p.u. voltage')
ax = gca;
ax.FontSize = 16;

subplot(4,1,2)
plot(vtrue.t*180/pi,'-o','linewidth',3);
hold on;
plot(vars.theta*180/pi,'-*','linewidth',3)
legend('true', 'linpf')
xlabel('bus number')
ylabel('bus angle [degree]')
ax = gca;
ax.FontSize = 16;

subplot(4,1,3)
plot(E*vtrue.t*180/pi,'-o','linewidth',3);
hold on;
plot(E*vars.theta*180/pi,'-*','linewidth',3)
legend('true', 'linpf')
xlabel('branch number')
ylabel('$\theta_{ft(\ell)}$ [degree]','Interpreter', 'Latex')
ax = gca;
ax.FontSize = 16;

subplot(4,1,4)
stem(find(bidx.pq),abs(residual(bidx.pq)))
hold on;
stem(find(bidx.pv),abs(residual(bidx.pv)))
stem(find(bidx.ref),abs(residual(bidx.ref)))
ylabel('|diag(v*)Yv - S*|')
xlabel('bus number')
legend('PQ', 'PV', 'REF')
ax = gca;
ax.FontSize = 16;