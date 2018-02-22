function S = powervectors(mpc,gmap)
% extract power vecotrs from mpc structure

S = struct();
define_constants;
S.Pg = gmap*mpc.gen(:,PG)/mpc.baseMVA; 
S.Pd = mpc.bus(:,PD)/mpc.baseMVA;
S.Qg = gmap*mpc.gen(:,QG)/mpc.baseMVA;
S.Qd = mpc.bus(:,QD)/mpc.baseMVA;