function vt_comp_plots(vars, vtrue, E, residual,bidx)
%%% assumption is that theta is given in radians!!!!
%% PQ
figure;
subplot(4,2,1)
plot(vtrue.v(bidx.pq),'-o','linewidth',3);
hold on;
if isreal(vars.v)
    plot(vars.v(bidx.pq),'-*','linewidth',3)
else
    plot(abs(vars.v(bidx.pq)),'-*','linewidth',3)
end
legend('true', 'linpf')
title('PQ buses')
xlabel('PQ bus number')
ylabel('p.u. voltage')
ax = gca;
ax.FontSize = 16;

subplot(4,2,3)
plot(vtrue.t(bidx.pq)*180/pi,'-o','linewidth',3);
hold on;
plot(vars.theta(bidx.pq)*180/pi,'-*','linewidth',3)
legend('true', 'linpf')
xlabel('PQ bus number')
ylabel('PQ bus angle [degree]')
ax = gca;
ax.FontSize = 16;

subplot(4,2,[5 6])
plot(E*vtrue.t*180/pi,'-o','linewidth',3);
hold on;
plot(E*vars.theta*180/pi,'-*','linewidth',3)
legend('true', 'linpf')
xlabel('branch number')
ylabel('$\theta_{ft(\ell)}$ [degree]','Interpreter', 'Latex')
ax = gca;
ax.FontSize = 16;

subplot(4,2,[7 8])
stem(find(bidx.pq),abs(residual(bidx.pq)))
hold on;
stem(find(bidx.pv),abs(residual(bidx.pv)))
stem(find(bidx.ref),abs(residual(bidx.ref)))
ylabel('|diag(v*)Yv - S*|')
xlabel('bus number')
legend('PQ', 'PV', 'REF')
ax = gca;
ax.FontSize = 16;

%% PV
subplot(4,2,2)
plot(imag(vtrue.sg(bidx.pv)),'-o','linewidth',3);
hold on;
plot(vars.Qpv,'-*','linewidth',3)
legend('true', 'linpf')
title('PV buses')
xlabel('PV bus number')
ylabel('Q_g [p.u.]')
ax = gca;
ax.FontSize = 16;

subplot(4,2,4)
plot(vtrue.t(bidx.pv)*180/pi,'-o','linewidth',3);
hold on;
plot(vars.theta(bidx.pv)*180/pi,'-*','linewidth',3)
legend('true', 'linpf')
xlabel('PV bus number')
ylabel('PV bus angle [degree]')
ax = gca;
ax.FontSize = 16;