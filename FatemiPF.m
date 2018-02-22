% linearized powerflow form Fatemi
% clear variables; close all;
define_constants;

lambda1 = 0.95;
lambda2 = 0.95;
casename = 'case118';
mpc = loadcase(casename);

N = size(mpc.bus,1);
M = size(mpc.branch,1);
G = size(mpc.gen,1);

% mpc.branch(:,TAP) = 0;
% mpc.branch(:,BR_B) = 0;
% mpc.branch(:,SHIFT) = 0;
% mpc.bus(:,GS) = 0;
% mpc.bus(:,BS) = 0;

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

if ~isempty(pv2pq) || ~isempty(ref2pq)
    fprintf('WARING: changed bus type')
end

Npq = sum(pq_idx);
Npv = sum(pv_idx);
Nref= sum(ref_idx);

[Y,Yf,Yt] = makeYbus(ext2int(mpc));
mpctmp = mpc;
mpctmp.bus = [mpc.bus(mpc.bus(:,BUS_TYPE) == REF,:); mpc.bus(mpc.bus(:,BUS_TYPE) == PV,:); mpc.bus(mpc.bus(:,BUS_TYPE) == PQ,:)];
Y2 = makeYbus(ext2int(mpctmp));

Gs = lambda1*real(Y2);
% Gs = Gs - sparse(1:N,1:N,Gs*ones(N,1),N,N);
Bs = lambda2*imag(Y2);
% Bs = Bs - sparse(1:N,1:N,Bs*ones(N,1),N,N);

Pg = gmap*mpc.gen(:,PG)/mpc.baseMVA; 
Pd = mpc.bus(:,PD)/mpc.baseMVA;
Qg = gmap*mpc.gen(:,QG)/mpc.baseMVA;
Qd = mpc.bus(:,QD)/mpc.baseMVA;

P  = -Pd;
P(pq_idx | pv_idx) = P(pq_idx | pv_idx) + Pg(pq_idx | pv_idx);
Q  = -Qd;
Q(pq_idx) = Q(pq_idx) + Qg(pq_idx);

v2pv  = vg(pv_idx).^2;
v2ref = vg(ref_idx).^2;
thetap_ref = v2ref.*mpc.bus(ref_idx,VA)*pi/180;

% structured matrices
mrange = 1:(Nref+Npv);
nrange = (Nref+Npv)+1:N;
v2m = [v2ref; v2pv];
Qn  = Q(pq_idx);
Ps  = [P(ref_idx);P(pv_idx);P(pq_idx)];
Qs  = [Q(ref_idx);Q(pv_idx);Q(pq_idx)];
% Bs = [B(ref_idx,:); B(~ref_idx,:)];
% Bs = [Bs(:,ref_idx), Bs(:,~ref_idx)];
% Bs = [B(ref_idx,:); B(pv_idx,:); B(pq_idx,:)];
% Gs = [G(ref_idx,:); G(pv_idx,:); G(pq_idx,:)];

Bnninv = Bs(nrange,nrange)\speye(Npq);
H11 = -(Bs(mrange,mrange) + Gs(mrange,nrange)*Bnninv*Gs(nrange,mrange));
H12 = -(Bs(mrange,nrange) + Gs(mrange,nrange)*Bnninv*Gs(nrange,nrange));
H21 = -(Bs(nrange,mrange) + Gs(nrange,nrange)*Bnninv*Gs(nrange,mrange));
H22 = -(Bs(nrange,nrange) + Gs(nrange,nrange)*Bnninv*Gs(nrange,nrange));

H = [H11, H12; H21, H22];

Lvm = Gs(mrange,mrange) - Gs(mrange,nrange)*Bnninv*Bs(nrange,mrange);
LQm = -Gs(mrange,nrange)*Bnninv;
Lvn = Gs(nrange,mrange) - Gs(nrange,nrange)*Bnninv*Bs(nrange,mrange);
LQn = -Gs(nrange,nrange)*Bnninv;

PvQ = [Lvm,LQm;Lvn,LQn]*[v2m;Qn];

thetap = [H11(2:end,2:end), H12(2:end,:); H21(:,2:end), H22]\(Ps(2:end) - PvQ(2:end));

v2n = -Bnninv*(Qn + Gs(nrange,mrange)*[0;thetap(1:Npv)] + Gs(nrange,nrange)*thetap(Npv+1:end) + Bs(nrange,mrange)*v2m);

theta = zeros(N,1);
theta(ref_idx) = thetap_ref./v2ref;
theta(pv_idx)  = theta(ref_idx) + thetap(1:Npv)./v2pv;
theta(pq_idx)  = theta(ref_idx) + thetap(Npv+1:end)./v2n;

% test_ans = Gs*[v2m;v2n] - Bs*[0;thetap] + Bs(:,1)*thetap_ref - P;
test_ans = Gs*[v2m;v2n] - Bs*[0;thetap] - Ps;

% Stest    = [Gs, -Bs; -Bs, -Gs]*[v2m;v2n;0;thetap];
% Ptest    = Stest(1:N);
% Qtest    = Stest(N+1:end);
% PvQtest  = [Lvm,LQm;Lvn,LQn]*[v2m;Qtest(nrange)];
% test_res = [H11(2:end,2:end), H12(2:end,:); H21(:,2:end), H22]\(Ptest(2:end) - PvQtest(2:end));
% figure;
% plot(test_res - thetap)
% return
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
plot(theta*180/pi,'-*','linewidth',3)
plot(mpcdc.bus(:,VA),'-x','linewidth',3);
legend('ac pf', 'Fatemi', 'dc pf')
xlabel('bus number')
ylabel('bus angle [degree]')
ax = gca;
ax.FontSize = 16;