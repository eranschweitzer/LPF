% check quality of linearization around correct solution
clear variables; close all;
define_constants;
casename = 'case118';
mpc = loadcase(casename);


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

[Yff,Yft,Ytf,Ytt,Ysh] = Yparts(mpc);
Gff = real(Yff); Bff=imag(Yff);
Gft = real(Yft); Bft=imag(Yft);
Gtf = real(Ytf); Btf=imag(Ytf);
Gtt = real(Ytt); Btt=imag(Ytt);

[BBUS, BF, PBUSINJ, PFINJ] = makeBdc(ext2int(mpc));


% solve AC
mpopt = mpoption;
mpopt.out.all = 0;
mpcac = runpf(mpc,mpopt);

%solved varialbes
theta = mpcac.bus(:,VA)*pi/180;
u     = log(mpcac.bus(:,VM));

breakpnt = mean(Gtt);
% weight = 1./(1 + (abs(Gtt)./breakpnt).^2);
weight = abs(1./(1 + 1i*2*0.7*abs(Gtt)./breakpnt + (1i*abs(Gtt)./breakpnt).^2));
% powerflow from linearized calculation
Pf  = exp(F*u).*(exp(F*u0).*(ones(M,1) + F*(u-u0)).*Gff.*weight + exp(T*u0).*Gft.*weight.*(ones(M,1) + T*(u-u0)) - exp(T*u0).*Bft.*(-E*theta));
Pf2 = exp(F*u).*exp(T*u0).*Bft.*(E*theta);
Pf3 = Bft.*(E*theta);
Sf = conj(exp(F*u).*(exp(F*u0).*(ones(M,1) + F*(u-u0)).*Yff + exp(T*u0).*Yft.*(ones(M,1) + T*(u-u0) + 1i*(-E*theta))));
St = conj(exp(T*u).*(exp(T*u0).*(ones(M,1) + T*(u-u0)).*Ytt + exp(F*u0).*Ytf.*(ones(M,1) + F*(u-u0) + 1i*(E*theta))));
DC = BF*theta + PFINJ;

% problem edge for case3375wp is 2068
% nf = mpc.branch(2068,1); nt = mpc.branch(2068,2); r = mpc.branch(2068,3); x = mpc.branch(2068,4); b = mpc.branch(2068,5);
% tneigh = find(E(:,nmap(nt))); fneigh = find(E(:,nmap(nf)));
parta = exp(F*u0).*(ones(M,1) + F*(u-u0));
partb = exp(T*u0).*(ones(M,1) + T*(u-u0));
partc = exp(T*u0).*(-1i*E*theta);

figure;
subplot(2,1,1)
plot(mpcac.branch(:,PF),'-o','linewidth',2);
hold on;
plot(Pf*mpc.baseMVA,'-*','linewidth',2);
plot(DC*mpc.baseMVA,'-*','linewidth',2);
xlabel('branch number')
ylabel('Pf')
yyaxis right;
plot(Pf*mpc.baseMVA - mpcac.branch(:,PF),'k--');
plot(DC*mpc.baseMVA - mpcac.branch(:,PF),'r--');
legend('ac pf', 'scaglione', 'dc', 'error', 'error dc')
ylabel('error')
ax = gca;
ax.FontSize = 16;

subplot(2,1,2)
plot(mpcac.branch(:,PT),'-o','linewidth',3);
hold on;
plot(real(St)*mpc.baseMVA,'-*','linewidth',3);
legend('ac pf', 'scaglione')
xlabel('branch number')
ylabel('Pt')
ax = gca;
ax.FontSize = 16;

%%
[er1,idx1]   = max(abs(Pf*mpc.baseMVA - mpcac.branch(:,PF)));
avg1         = mean(abs(Pf*mpc.baseMVA - mpcac.branch(:,PF)));
std1         = std(abs(Pf*mpc.baseMVA - mpcac.branch(:,PF)));
[er2,idx2]   = max(abs(Pf2*mpc.baseMVA - mpcac.branch(:,PF)));
[er3,idx3]   = max(abs(Pf3*mpc.baseMVA - mpcac.branch(:,PF)));
[erdc,idxdc] = max(abs(DC*mpc.baseMVA  - mpcac.branch(:,PF)));
fprintf('max(|error1|) = %0.2f MW (%0.2f%%), mean=%0.2f, std=%0.2f\n',er1,100*er1/abs(mpcac.branch(idx1,PF)),avg1,std1)
fprintf('max(|error dc|) = %0.2f MW (%0.2f%%)\n',erdc,100*erdc/abs(mpcac.branch(idxdc,PF))) 
fprintf('max(|error2|) = %0.2f MW (%0.2f%%)\n',er2,100*er2/abs(mpcac.branch(idx2,PF)))
fprintf('max(|error3|) = %0.2f MW (%0.2f%%)\n',er3,100*er3/abs(mpcac.branch(idx3,PF)))
%%
nf = mpc.branch(idx1,1); nt = mpc.branch(idx1,2); r = mpc.branch(idx1,3); x = mpc.branch(idx1,4); b = mpc.branch(idx1,5);
tneigh = find(E(:,nmap(nt))); fneigh = find(E(:,nmap(nf)));