%Coffrin model
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
E  = sparse([1:M,1:M].',[full(nmap(mpc.branch(:,F_BUS)));full(nmap(mpc.branch(:,T_BUS)))],[ones(M,1);-1*ones(M,1)],M,N);
gbus  = full(nmap(mpc.gen(:,GEN_BUS)));
gmap  = sparse(gbus,1:G,(mpc.gen(:,GEN_STATUS) > 0),N,G);
bus_with_ongen = sum(gmap,2) > 0;
ref  = find(mpc.bus(:,BUS_TYPE) == 3);
pq_idx  = mpc.bus(:,BUS_TYPE) == PQ;
pv_idx  = mpc.bus(:,BUS_TYPE) == PV;
ref_idx = mpc.bus(:,BUS_TYPE) == REF;


% target voltage. Everything is 1 except for the PV and ref buses
vg = ones(N,1);
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

h = 20;
anglim = pi/3;
d = 2*anglim/(h+1);

[Y,Yf,Yt] = makeYbus(ext2int(mpc));

gft = zeros(M,1);
bft = zeros(M,1);
gtf = zeros(M,1);
btf = zeros(M,1);
gff = zeros(M,1);
gtt = zeros(M,1);
bff = zeros(M,1);
btt = zeros(M,1);

for k = 1:M
    fbus = full(nmap(mpc.branch(k,F_BUS)));
    tbus = full(nmap(mpc.branch(k,T_BUS)));
    
    yft  = -Y(fbus,tbus);
    ytf  = -Y(tbus,fbus);
    
    yff  = Yf(k,fbus);% - 1i*mpc.branch(k,BR_B)/(2*tau^2);
    ytt  = Yt(k,tbus);% - 1i*mpc.branch(k,BR_B)/2;
    
    gft(k) = real(yft);
    bft(k) = imag(yft);
    
    gtf(k) = real(ytf);
    btf(k) = imag(ytf);
    
    gff(k) = real(yff);
    bff(k) = imag(yff);
    
    gtt(k) = real(ytt);
    btt(k) = imag(ytt);
    
end


%===============
% variables
%===============
theta = sdpvar(N,1);
phi   = sdpvar(N,1);
cosf  = sdpvar(M,1);
cost  = sdpvar(M,1);
pf    = sdpvar(M,1);
pt    = sdpvar(M,1);
qf    = sdpvar(M,1);
qt    = sdpvar(M,1);
qfd    = sdpvar(M,1);
qtd    = sdpvar(M,1);

vf = vg(full(nmap(mpc.branch(:,F_BUS))));
vt = vg(full(nmap(mpc.branch(:,T_BUS))));


Pg = gmap*mpc.gen(:,PG)/mpc.baseMVA; 
Pd = mpc.bus(:,PD)/mpc.baseMVA;
Qg = gmap*mpc.gen(:,QG)/mpc.baseMVA;
Qd = mpc.bus(:,QD)/mpc.baseMVA;

% P  = -(Pd + mpc.bus(:,GS)/mpc.baseMVA);
P  = -Pd;
P(pq_idx | pv_idx) = P(pq_idx | pv_idx) + Pg(pq_idx | pv_idx);
% Q  = -(Qd - mpc.bus(:,BS)/mpc.baseMVA);
Q  = -Qd;
Q(pq_idx) = Q(pq_idx) + Qg(pq_idx);

sumF = sparse(full(nmap(mpc.branch(:,F_BUS))),1:M,1,N,M);
sumT = sparse(full(nmap(mpc.branch(:,T_BUS))),1:M,1,N,M);

constraints = [...
pf  ==  vf.^2.*gft - vf.*vt.*(gft.*cosf + bft.*E*theta),...
pt  ==  vt.^2.*gtf - vt.*vf.*(gtf.*cost - btf.*E*theta),...
qf  == -vf.^2.*bft - vf.*vt.*(gft.*E*theta - bft.*cosf),...
qt  == -vt.^2.*btf - vt.*vf.*(-gtf.*E*theta - btf.*cost),...
qfd == -vf.*bft.*E*phi - E*vg.*bft.*phi(full(nmap(mpc.branch(:,F_BUS)))),...
qtd ==  vt.*btf.*E*phi + E*vg.*btf.*phi(full(nmap(mpc.branch(:,T_BUS)))),...
P(~ref_idx) == sumF(~ref_idx,:)*pf + sumT(~ref_idx,:)*pt,...
Q(pq_idx)   == sumF(pq_idx,:)*(qf + qfd) + sumT(pq_idx,:)*(qt + qtd),...
theta(ref_idx) == mpc.bus(ref_idx,VA)*pi/180,...
phi(~pq_idx)   == 0, phi >= -vg,...
cos(anglim) <= cosf <= 1, cos(anglim) <= cost <= 1,...
cost == cosf,...
];

for t = 1:h
    constraints = [constraints,...
        cosf <= -sin(t*d-anglim)*( E*theta - t*d + anglim) + cos(t*d -anglim)];%,...
%         cost <= -sin(t*d-anglim)*(-E*theta - t*d + anglim) + cos(t*d - anglim)];
end

objective = -sum(cosf);

res = optimize(constraints,objective);

%%
mpopt = mpoption;
mpopt.out.all = 0;
mpcac = runpf(mpc,mpopt);
mpcac.If = Yf*(mpcac.bus(:,VM).*exp(1i*pi/180*mpcac.bus(:,VA)));
mpcac.It = Yt*(mpcac.bus(:,VM).*exp(1i*pi/180*mpcac.bus(:,VA)));
mpcac.Pf = mpcac.bus(nmap(mpcac.branch(:,1)),VM).*exp(1i*pi/180*mpcac.bus(nmap(mpcac.branch(:,1)),VA)).*conj(mpcac.If);
mpcdc = rundcpf(mpc,mpopt);

%%
ang.abserr  = abs(mpcac.bus(:,VA)*pi/180 - value(theta));
ang.corr = corr(mpcac.bus(:,VA)*pi/180,value(theta));
ang.muerr   = mean(ang.abserr);
ang.maxerr  = max(ang.abserr);

p.f.abserr    = abs(mpcac.branch(:,PF) - value(pf)*100);
p.f.corr      = corr(mpcac.branch(:,PF),value(pf)*100);
p.f.muerr     = mean(p.f.abserr);
p.f.maxerr    = max(p.f.abserr);

p.t.abserr    = abs(mpcac.branch(:,PT) - value(pt)*100);
p.t.corr      = corr(mpcac.branch(:,PT),value(pt)*100);
p.t.muerr     = mean(p.t.abserr);
p.t.maxerr    = max(p.t.abserr);


%%
figure;
plot(mpcac.bus(:,VA),'-o','linewidth',3);
hold on;
plot(value(theta)*180/pi,'-*','linewidth',3)
plot(mpcdc.bus(:,VA),'-x','linewidth',3);
legend('ac pf', 'coffrin', 'dc pf')
xlabel('bus number')
ylabel('bus angle [degree]')
ax = gca;
ax.FontSize = 16;