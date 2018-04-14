function vtrue = vtrue_struct(mpc, opmats, Sb, Sp)
define_constants;
v2struct(opmats);
%% AC solution
mpopt   = mpoption; mpopt.out.all = 0;
mpcac   = runpf(mpc,mpopt);
vtrue   = struct();
vtrue.v = mpcac.bus(:,VM);        %true voltage magnitude solution
vtrue.t = mpcac.bus(:,VA)*pi/180; %true voltage angle solution
vtrue.phi = 0.5*(E*vtrue.t).^2;       %true phi (1/2)*(delta(theta))^2
vtrue.pf  = mpcac.branch(:,PF);     %true real power at from bus [MW]
vtrue.pt  = mpcac.branch(:,PT);     %true real power at to bus [MW]
vtrue.qf  = mpcac.branch(:,QF);     %true reactive power at from bus [MVAr]
vtrue.qt  = mpcac.branch(:,QT);     %true reactive power at to bus [MVAr]
vtrue.sg  = gmap*(mpcac.gen(:,PG) + 1i*mpcac.gen(:,QG))/baseMVA;
vtrue.residual = pfresidual(vtrue.v.*exp(1i*vtrue.t), myMakeYbus(F,T,Sb), real(vtrue.sg)-Sp.Pd, imag(vtrue.sg)-Sp.Qd );