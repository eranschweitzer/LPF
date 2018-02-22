clear variables; close all;
define_constants;
casename = 'case_ACTIVSg2000';
mpc = loadcase(casename);

% mpc.bus(:,BS) = 0;
% mpc.bus(:,GS) = 0;
% mpc.branch(:,BR_B) = 0;

N = size(mpc.bus,1);
M = size(mpc.branch,1);
G = size(mpc.gen,1);


nmap  = sparse(mpc.bus(:,BUS_I),1,1:N);
E  = sparse([1:M,1:M].',[full(nmap(mpc.branch(:,F_BUS)));full(nmap(mpc.branch(:,T_BUS)))],[ones(M,1);-1*ones(M,1)],M,N);
F  = sparse(1:M,full(nmap(mpc.branch(:,F_BUS))),1,M,N);
T  = sparse(1:M,full(nmap(mpc.branch(:,T_BUS))),1,M,N);
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

upv = log(vg(pv_idx));
theta_ref = mpc.bus(ref_idx,VA)*pi/180;
uref = log(vg(ref_idx));
u0 = uref*ones(N,1);
u0(pv_idx)  = upv;
u0(ref_idx) = uref;

% load and generation
Pg = gmap*mpc.gen(:,PG)/mpc.baseMVA; 
Pd = mpc.bus(:,PD)/mpc.baseMVA;
Qg = gmap*mpc.gen(:,QG)/mpc.baseMVA;
Qd = mpc.bus(:,QD)/mpc.baseMVA;

P  = Pg-Pd;
% P(pq_idx | pv_idx) = P(pq_idx | pv_idx) + Pg(pq_idx | pv_idx);

Q  = -Qd;
Q(pq_idx) = Q(pq_idx) + Qg(pq_idx);

% admitances
[Yff,Yft,Ytf,Ytt,Ysh] = Yparts(mpc);
% Gff = real(Yff); Bff=imag(Yff);
% Gft = real(Yft); Bft=imag(Yft);
% Gtf = real(Ytf); Btf=imag(Ytf);
% Gtt = real(Ytt); Btt=imag(Ytt);

breakpnt_G = mean(abs(real(Ytt)));
wG = abs(1./(1 + 1i*2*0.7*abs(real(Ytt))./breakpnt_G + (1i*abs(real(Ytt))./breakpnt_G).^2));

breakpnt_B = mean(abs(imag(Ytt)))+2*std(abs(imag(Ytt)));
wB = abs(1./(1 + 1i*2*0.7*abs(imag(Ytt))./breakpnt_B + (1i*abs(imag(Ytt))./breakpnt_B).^2));

Rxft = -real(Yft)./imag(Yft);
Rxtf = -real(Ytf)./imag(Ytf);
% form sub-matrices
delta2 = -2.5e-3;%0.5*(sum(Pd)*median(imag(Ytf))/M^2).^2;
Atr = ( F'*sparse(1:M,1:M,imag(Yft).*exp(T*u0),M,M) - T'*sparse(1:M,1:M,imag(Ytf).*exp(F*u0),M,M))*E;
Atm = (-F'*sparse(1:M,1:M,real(Yft).*exp(T*u0),M,M) + T'*sparse(1:M,1:M,real(Ytf).*exp(F*u0),M,M))*E;
Aur = F'*(sparse(1:M,1:M, wG.*real(Yff).*exp(F*u0),M,M)*F + sparse(1:M,1:M, wG.*real(Yft).*exp(T*u0),M,M)*T) + ...
      T'*(sparse(1:M,1:M, wG.*real(Ytt).*exp(T*u0),M,M)*T + sparse(1:M,1:M, wG.*real(Ytf).*exp(F*u0),M,M)*F) + ...
      sparse(1:N,1:N,real(Ysh).*exp(u0) + P.*exp(-u0),N,N);
Aum = F'*(sparse(1:M,1:M, wB.*imag(Yff).*exp(F*u0),M,M)*F + sparse(1:M,1:M, wB.*imag(Yft).*exp(T*u0),M,M)*T) + ...
      T'*(sparse(1:M,1:M, wB.*imag(Ytt).*exp(T*u0),M,M)*T + sparse(1:M,1:M, wB.*imag(Ytf).*exp(F*u0),M,M)*F) + ...
      sparse(1:N,1:N,imag(Ysh).*exp(u0) - Q.*exp(-u0),N,N);
% br   =  P.*exp(-2*u0) - real(Ysh) - ...
%     (F'*(real(Yff) + real(Yft).*exp(-E*u0)) + T'*(real(Ytt) + real(Ytf).*exp(E*u0)));
% bm   = -Q.*exp(-2*u0) - imag(Ysh) - ...
%     (F'*(imag(Yff) + imag(Yft).*exp(-E*u0)) + T'*(imag(Ytt) + imag(Ytf).*exp(E*u0)));

br =  P.*exp(-u0) - real(Ysh).*exp(u0) - ...
      F'*(wG.*real(Yff).*exp(F*u0) + wG.*real(Yft).*exp(T*u0).*(1 + delta2)) - ...
      T'*(wG.*real(Ytt).*exp(T*u0) + wG.*real(Ytf).*exp(F*u0).*(1 + delta2));
bm = -Q.*exp(-u0) - imag(Ysh).*exp(u0) - ...
      F'*(wB.*imag(Yff).*exp(F*u0) + wB.*imag(Yft).*exp(T*u0).*(1 + delta2)) - ...
      T'*(wB.*imag(Ytt).*exp(T*u0) + wB.*imag(Ytf).*exp(F*u0).*(1 - delta2));

Ipv  = sparse(1:Npv,find(pv_idx),1,Npv,N);
Upv  = sparse(find(pv_idx),1:Npv, exp(-u0(pv_idx)), N, Npv);

Iref  = sparse(1:Nref,find(ref_idx),1, Nref, N);
Uref  = sparse(find(ref_idx),1:Nref, exp(-u0(ref_idx)), N, Nref);

A = [Atr,Aur,sparse(N,Npv),-Uref,sparse(N,Nref);
     Atm,Aum,Upv,sparse(N,Nref),Uref;
     sparse(Npv,N),Ipv,sparse(Npv,Npv+2*Nref);
     Iref,sparse(Nref,N+Npv+2*Nref);
     sparse(Nref,N),Iref,sparse(Nref,Npv+2*Nref);
     ];
b = [br;bm;zeros(Npv,1);theta_ref;zeros(Nref,1)];

% solve
x = A\b;

% theta = x(1:N);
% delta2 = 0.5*(E*theta).^2;
% br =  P.*exp(-u0) - real(Ysh).*exp(u0) - ...
%       F'*(wG.*real(Yff).*exp(F*u0) + wG.*real(Yft).*exp(T*u0).*(1 + delta2)) - ...
%       T'*(wG.*real(Ytt).*exp(T*u0) + wG.*real(Ytf).*exp(F*u0).*(1 + delta2));
% bm = -Q.*exp(-u0) - imag(Ysh).*exp(u0) - ...
%       F'*(wB.*imag(Yff).*exp(F*u0) + wB.*imag(Yft).*exp(T*u0).*(1 + delta2)) - ...
%       T'*(wB.*imag(Ytt).*exp(T*u0) + wB.*imag(Ytf).*exp(F*u0).*(1 - delta2));
% b = [br;bm;zeros(Npv,1);theta_ref;zeros(Nref,1)];
% x = A\b;

%parse result
theta = x(1:N);
u     = x(N+1:2*N) + u0;
v     = exp(u);
Qpv   = x(2*N+1:2*N+Npv);
Pref  = x(2*N+Npv+1:2*N+Npv+Nref);
Qref  = x(end-Nref+1:end);

Pf = exp(F*u).*(exp(F*u0).*(ones(M,1) + F*(u-u0)).*wG.*real(Yff) + ...
                exp(T*u0).*(ones(M,1) + T*(u-u0) - delta2).*wG.*real(Yft) + ...
               -exp(T*u0).*imag(Yft).*(-E*theta));

Pt = exp(T*u).*(exp(T*u0).*(ones(M,1) + T*(u-u0)).*wG.*real(Ytt) + ...
                exp(F*u0).*(ones(M,1) + F*(u-u0) - delta2).*wG.*real(Ytf) + ...
               -exp(F*u0).*imag(Ytf).*(E*theta));
Sf = conj(exp(F*u).*(exp(F*u0).*(ones(M,1) + F*(u-u0)).*Yff + exp(T*u0).*Yft.*(ones(M,1) + T*(u-u0) + 1i*(-E*theta))));
St = conj(exp(T*u).*(exp(T*u0).*(ones(M,1) + T*(u-u0)).*Ytt + exp(F*u0).*Ytf.*(ones(M,1) + F*(u-u0) + 1i*(E*theta))));
%%
mpopt = mpoption;
mpopt.out.all = 0;
mpcac = runpf(mpc,mpopt);
mpcdc = rundcpf(mpc,mpopt);
%% error
scag.Pf = Pf*mpc.baseMVA - mpcac.branch(:,PF);
scag.Pt = Pt*mpc.baseMVA - mpcac.branch(:,PT);
dc.Pf   = mpcdc.branch(:,PF) - mpcac.branch(:,PF);
dc.Pt   = mpcdc.branch(:,PT) - mpcac.branch(:,PT);

[er1,idx1]   = max(abs(scag.Pf));
[er2,idx2]   = max(abs(dc.Pf));
fprintf('max(|scaglione|) = %0.2f MW (%0.2f%%)\n',er1,100*er1/abs(mpcac.branch(idx1,PF)))
fprintf('max(|dc|) = %0.2f MW (%0.2f%%)\n',er2,100*er2/abs(mpcac.branch(idx2,PF)))
%%
nf = mpc.branch(idx1,1); nt = mpc.branch(idx1,2); r = mpc.branch(idx1,3); x = mpc.branch(idx1,4); b = mpc.branch(idx1,5);
tneigh = find(E(:,nmap(nt))); fneigh = find(E(:,nmap(nf)));
%%
figure;
plot(mpcac.bus(:,VM),'-o','linewidth',3);
hold on;
plot(v,'-*','linewidth',3)
legend('ac pf', 'scaglione')
xlabel('bus number')
ylabel('p.u. voltage')
ax = gca;
ax.FontSize = 16;

figure;
plot(mpcac.bus(:,VA),'-o','linewidth',3);
hold on;
plot(theta*180/pi,'-*','linewidth',3)
plot(mpcdc.bus(:,VA),'-x','linewidth',3);
legend('ac pf', 'scaglione', 'dc pf')
xlabel('bus number')
ylabel('bus angle [degree]')
ax = gca;
ax.FontSize = 16;

figure;
subplot(2,1,1)
plot(mpcac.branch(:,PF),'-o','linewidth',3);
hold on;
plot(Pf*mpc.baseMVA,'-*','linewidth',3);
plot(mpcdc.branch(:,PF),'-x','linewidth',3);
xlabel('branch number')
ylabel('Pf')
yyaxis right;
plot(scag.Pf,'k--');
hold on;
plot(mpcdc.branch(:,PF) - mpcac.branch(:,PF),'r--');
legend('ac pf', 'scaglione', 'dc pf','scaglione error', 'dc pf error')
ylabel('error')
ax = gca;
ax.FontSize = 16;

subplot(2,1,2)
plot(mpcac.branch(:,PT),'-o','linewidth',3);
hold on;
plot(Pt*mpc.baseMVA,'-*','linewidth',3);
plot(mpcdc.branch(:,PT),'-x','linewidth',3);
legend('ac pf', 'scaglione', 'dc pf')
xlabel('branch number')
ylabel('Pt')
ax = gca;
ax.FontSize = 16;