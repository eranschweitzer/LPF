%LPPF
clear variables; close all;
define_constants;
casename = 'case118';
mpc = loadcase(casename);

% mpc.branch(:,BR_B) = 0;
% mpc.branch(:,TAP) = 0;
% mpc.branch(:,SHIFT) = 0;
% mpc.bus(:,BS) = 0;
% mpc.bus(:,GS) = 0;

N = size(mpc.bus,1);
M = size(mpc.branch,1);
G = size(mpc.gen,1);


nmap  = sparse(mpc.bus(:,BUS_I),1,1:N);
gbus  = full(nmap(mpc.gen(:,GEN_BUS)));
gmap  = sparse(gbus,1:G,(mpc.gen(:,GEN_STATUS) > 0),N,G);
bus_with_ongen = sum(gmap,2) > 0;
ref  = find(mpc.bus(:,BUS_TYPE) == 3);
pq_idx  = mpc.bus(:,BUS_TYPE) == PQ;
pv_idx  = mpc.bus(:,BUS_TYPE) == PV;
ref_idx = mpc.bus(:,BUS_TYPE) == REF;


vg = zeros(N,1);
vg(gbus) = mpc.gen(:,VG); 

% if a bus is specified as PV but has no generator attached change it to PQ
% same for ref bus (though this really shouldn't happen normally)
pv2pq  = find(pv_idx  & ~bus_with_ongen);
pv_idx(pv2pq)   = false;
pq_idx(pv2pq)   = true;

ref2pq = find(ref_idx & ~bus_with_ongen);
ref_idx(ref2pq) = false;
ref_idx(ref2pq) = true;

Npq = sum(pq_idx);
Npv = sum(pv_idx);
Nref= sum(ref_idx);

[Y,Yf,Yt] = makeYbus(ext2int(mpc));
ysh = Y*ones(N,1);
Ytr = Y - sparse(1:N,1:N,ysh,N,N);

Gtr = real(Ytr);
Btr = imag(Ytr);
gsh = real(ysh);
bsh = imag(ysh);

Pg = gmap*mpc.gen(:,PG)/mpc.baseMVA; 
Pd = mpc.bus(:,PD)/mpc.baseMVA;
Qg = gmap*mpc.gen(:,QG)/mpc.baseMVA;
Qd = mpc.bus(:,QD)/mpc.baseMVA;

upv = log(vg(pv_idx));
theta_ref = mpc.bus(ref_idx,VA)*pi/180;
uref = log(vg(ref_idx));

P  = -Pd;
P(pq_idx | pv_idx) = P(pq_idx | pv_idx) + Pg(pq_idx | pv_idx);
Q  = -Qd;
Q(pq_idx) = Q(pq_idx) + Qg(pq_idx);

dthetalim = 60*pi/180;
htheta    = 11;
theta_step = 2*dthetalim/(htheta+1);
umin = -0.3;
umax = 0.3;
hu   = 7;
ustep = (umax - umin)/(hu+1);
slope_u = (exp(-2*umax) - exp(-2*umin))/(umax - umin);

u = sdpvar(N,1);
theta = sdpvar(N,1);
w = sdpvar(N,1);
z = sdpvar(M,1);
pref = sdpvar(Nref,1);
qref = sdpvar(Nref,1);
qpv  = sdpvar(Npv,1);

zi = zeros(2*M,1);
zj = zeros(2*M,1);
gzv = zeros(2*M,1);
bzv = zeros(2*M,1);
ptr = 1;
for k = 1:M
    fbus = full(nmap(mpc.branch(k,F_BUS)));
    tbus = full(nmap(mpc.branch(k,T_BUS)));
    zi(ptr) = fbus;
    zj(ptr) = k;
    gzv(ptr)= Gtr(fbus,tbus);
    bzv(ptr)= Btr(fbus,tbus);
    ptr = ptr + 1;
    zi(ptr) = tbus;
    zj(ptr) = k;
    gzv(ptr)= -Gtr(tbus,fbus);
    bzv(ptr)= -Btr(tbus,fbus);
    ptr = ptr + 1;
end
Gz = sparse(zi,zj,gzv,N,M);
Bz = sparse(zi,zj,bzv,N,M);
E  = sparse([1:M,1:M].',[full(nmap(mpc.branch(:,F_BUS)));full(nmap(mpc.branch(:,T_BUS)))],[ones(M,1);-1*ones(M,1)],M,N);

constraints = [...
    -P(pq_idx).*w(pq_idx) + Gtr(pq_idx,:)*u - Btr(pq_idx,:)*theta - Gz(pq_idx,:)*z == -gsh(pq_idx),...
    Gtr(pv_idx,:)*u - Btr(pv_idx,:)*theta - Gz(pv_idx,:)*z == P(pv_idx).*exp(-2*upv) - gsh(pv_idx),...
    -pref.*exp(-2*uref) + Gtr(ref_idx,:)*u - Btr(ref_idx,:)*theta - Gz(ref_idx,:)*z == -Pd(ref_idx).*exp(-2*uref) - gsh(ref_idx),...
    Q(pq_idx,:).*w(pq_idx) + Btr(pq_idx,:)*u + Gtr(pq_idx,:)*theta - Bz(pq_idx,:)*z == -bsh(pq_idx),...
    exp(-2*upv).*qpv + Btr(pv_idx,:)*u + Gtr(pv_idx,:)*theta - Bz(pv_idx,:)*z == Qd(pv_idx).*exp(-2*upv) - bsh(pv_idx),...
    exp(-2*uref).*qref + Btr(ref_idx,:)*u + Gtr(ref_idx,:)*theta - Bz(ref_idx,:)*z == Qd(ref_idx).*exp(-2*uref) - bsh(ref_idx),...
    u(pv_idx) == upv,...
    u(ref_idx) == uref,...
    theta(ref_idx) == theta_ref,...
    w <= slope_u*u - slope_u*umin + exp(-2*umin),...
%     w == 1,...
    z <= dthetalim^2/2,...
    ];

for t = 1:htheta
    constraints = [constraints, ...
        z >= (-dthetalim + t*theta_step)^2/2 + (-dthetalim + t*theta_step)*(E*theta - (-dthetalim + t*theta_step))];
end

for t = 1:hu
    constraints = [constraints,...
        w >= -2*exp(-2*(umin + t*ustep))*(u - (umin + t*ustep)) + exp(-2*(umin+t*ustep))];
end

objective = sum(w) + sum(z);

diagnostics = optimize(constraints,objective);
%%
mpopt = mpoption;
mpopt.out.all = 0;
mpcac = runpf(mpc,mpopt);
mpcac.If = Yf*(mpcac.bus(:,VM).*exp(1i*pi/180*mpcac.bus(:,VA)));
mpcac.Pf = mpcac.bus(nmap(mpcac.branch(:,1)),VM).*exp(1i*pi/180*mpcac.bus(nmap(mpcac.branch(:,1)),VA)).*conj(mpcac.If);
mpcdc = rundcpf(mpc,mpopt);

%%
figure;
plot(mpcac.bus(:,VA),'-o','linewidth',3);
hold on;
plot(value(theta)*180/pi,'-*','linewidth',3)
plot(mpcdc.bus(:,VA),'-x','linewidth',3);
legend('ac pf', 'scaglione', 'dc pf')
xlabel('bus number')
ylabel('bus angle [degree]')
ax = gca;
ax.FontSize = 16;

fprintf('sum of absolute error (degree): %0.3f\n', sum(abs(mpcac.bus(:,VA) - value(theta)*180/pi)))