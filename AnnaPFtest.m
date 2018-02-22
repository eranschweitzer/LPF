clear variables; close all;
define_constants;
ufp = 1;%0.75;
ufq = 1;%0.25;
ufs = 1;%0.75;%0.5;
lambda1 = 0.95;
lambda2 = 0.95;
casename = '../cycles/ERCOT_P16SP.m';
casename = 'case_ACTIVSg2000';
version = 2;
mpc = loadcase(casename);

% mpc.bus(:,BS) = 0;
% mpc.bus(:,GS) = 0;
% mpc.branch(:,BR_B) = 0;

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

% ytr   = (mpc.branch(:,BR_R) + 1i*mpc.branch(:,BR_X)).^(-1);
% yc    = 1i*mpc.branch(:,BR_B);
% ysh   = (mpc.branch(:,GS) + 1i*mpc.branch(:,BS))./mpc.baseMVA;
% tap   = mpc.branch(:,TAP);
% shift = mpc.branch(:,SHIFT);

[Y,Yf,Yt] = makeYbus(ext2int(mpc));
ysh = Y*ones(N,1);
Ytr = Y - sparse(1:N,1:N,ysh,N,N);
% Ytr = Y - sparse(1:N,1:N,diag(Y),N,N);
% Ytr = Ytr - sparse(1:N,1:N,Ytr*ones(N,1),N,N);
% ysh = diag(Y - Ytr);

Gtr = real(Ytr);
Btr = imag(Ytr);
gsh = real(ysh);
bsh = imag(ysh);
% ufs = ufs./(max(abs(bsh),1));

Pg = gmap*mpc.gen(:,PG)/mpc.baseMVA; 
Pd = mpc.bus(:,PD)/mpc.baseMVA;
Qg = gmap*mpc.gen(:,QG)/mpc.baseMVA;
Qd = mpc.bus(:,QD)/mpc.baseMVA;

if ~strcmp(version,'fatemi')
    upv = log(vg(pv_idx));
    % upv = log(mpc.bus(pv_idx,VM));
    theta_ref = mpc.bus(ref_idx,VA)*pi/180;
    uref = log(vg(ref_idx));
    % uref = log(mpc.bus(ref_idx,VM));
else
    v2pv  = vg(pv_idx).^2;
    v2ref = vg(ref_idx).^2;
    theta_ref = mpc.bus(ref_idx,VA)*pi/180;
end

if all(version == 1) || strcmp(version,'fatemi')
    P  = -Pd;
    P(pq_idx | pv_idx) = P(pq_idx | pv_idx) + Pg(pq_idx | pv_idx);
    Q  = -Qd;
    Q(pq_idx) = Q(pq_idx) + Qg(pq_idx);
elseif all(version == 2)
    P  = zeros(N,1);
    P(pq_idx)  = Pg(pq_idx) - Pd(pq_idx);
    P(pv_idx)  = exp(-upv).*(Pg(pv_idx) - Pd(pv_idx));
    P(ref_idx) = -exp(-uref).*Pd(ref_idx);
    Q = zeros(N,1);
    Q(pq_idx) = Qg(pq_idx) - Qd(pq_idx);
    Q(pv_idx) = -exp(-upv).*Qd(pv_idx);
    Q(ref_idx) = -exp(-uref).*Qd(ref_idx);
elseif all(version == 3)
    P  = zeros(N,1);
    P(pq_idx)  = Pg(pq_idx) - Pd(pq_idx);
    P(pv_idx)  = exp(-2*upv).*(Pg(pv_idx) - Pd(pv_idx));
    P(ref_idx) = -exp(-2*uref).*Pd(ref_idx);
    Q = zeros(N,1);
    Q(pq_idx) = Qg(pq_idx) - Qd(pq_idx);
    Q(pv_idx) = -exp(-2*upv).*Qd(pv_idx);
    Q(ref_idx) = -exp(-2*uref).*Qd(ref_idx);
end



Ptilde = P;
Ptilde(pv_idx | ref_idx) = 0;
Qtilde = Q;
Qtilde(pv_idx | ref_idx) = 0;

Ipv  = sparse(1:Npv,find(pv_idx),1,Npv,N);
if all(version == 1)
    Upv  = sparse(find(pv_idx),1:Npv, ufq*upv - 1, N, Npv);
elseif all(version == 2)
    Upv  = sparse(find(pv_idx),1:Npv, -exp(-upv), N, Npv);
elseif all(version == 3)
    Upv  = sparse(find(pv_idx),1:Npv, -exp(-2*upv), N, Npv);
end
% Upv  = sparse(find(pv_idx),1:sum(pv_idx),  -1, N, sum(pv_idx));

Iref = sparse(1:Nref,find(ref_idx),1,Nref,N);
if all(version == 1)
    Urefp  = sparse(find(ref_idx),1:Nref, ufp*uref - 1, N, Nref);
    Urefq  = sparse(find(ref_idx),1:Nref, ufq*uref - 1, N, Nref);
elseif all(version == 2)
    Urefp  = sparse(find(ref_idx),1:Nref, -exp(-uref), N, Nref);
    Urefq  = sparse(find(ref_idx),1:Nref, -exp(-uref), N, Nref);
elseif all(version == 3)
    Urefp  = sparse(find(ref_idx),1:Nref, -exp(-2*uref), N, Nref);
    Urefq  = sparse(find(ref_idx),1:Nref, -exp(-2*uref), N, Nref);
end

% Uref  = sparse(find(ref_idx),1:sum(ref_idx), -1, N, sum(ref_idx));

% remove slack bus
% P(ref) = [];
% Q(ref) = [];
% Gtr(:,ref) = []; Gtr(ref,:) = [];
% Btr(:,ref) = []; Btr(ref,:) = [];
% gsh(ref) = [];
% bsh(ref) = [];
% 
% A = 2*sparse(1:2*(N-1),[1:N-1,1:N-1], [P;Q], 2*(N-1),2*(N-1)) + [-Gtr, Btr; Btr, Gtr];
% b = [P;Q] - [gsh;-bsh];
% x = A\b;
% u = x(1:N-1);
% v = exp(u);
% theta = x(N:end);
% return

if all(version == 1)
    Puterm =  sparse(1:N,1:N,ufp*P+ufs.*gsh);
    Quterm =  sparse(1:N,1:N,ufq*Q-ufs.*bsh);
elseif all(version == 2)
    Puterm =  sparse(1:N,1:N,ufp*Ptilde+ufs.*gsh);
    Quterm =  sparse(1:N,1:N,ufq*Qtilde-ufs.*bsh);
elseif all(version == 3)
    Puterm =  sparse(1:N,1:N,ufp*Ptilde+ufs.*gsh.*pq_idx);
    Quterm =  sparse(1:N,1:N,ufq*Qtilde-ufs.*bsh.*pq_idx);
end

if ~strcmp(version,'fatemi')
    A = [-Btr, Gtr+Puterm, sparse(N,Npv), Urefp, sparse(N,Nref);
         -Gtr, -Btr+Quterm, Upv, sparse(N,Nref), Urefq;
         sparse(Npv,N), Ipv, sparse(Npv,Npv+2*Nref);
         Iref, sparse(Nref, N+Npv+2*Nref);
         sparse(Nref,N), Iref, sparse(Nref,Npv+2*Nref)];

    b = [P-ufs.*gsh;Q+ufs.*bsh;upv;theta_ref;uref];
else
    A2 = [-lambda1*Btr, lambda1*Gtr+sparse(1:N,1:N,lambda1*gsh,N,N), sparse(N,Npv), -Iref', sparse(N,Nref);
        -lambda2*Gtr, -lambda2*Btr+sparse(1:N,1:N,-lambda2*bsh,N,N), -Ipv', sparse(N,Nref), -Iref';
        sparse(Npv,N), Ipv, sparse(Npv,Npv+2*Nref);
         Iref, sparse(Nref, N+Npv+2*Nref);
         sparse(Nref,N), Iref, sparse(Nref,Npv+2*Nref)];
    A = [-lambda1*imag(Y), lambda1*real(Y), sparse(N,Npv), -Iref', sparse(N,Nref);
        -lambda2*real(Y), -lambda2*imag(Y), -Ipv', sparse(N,Nref), -Iref';
        sparse(Npv,N), Ipv, sparse(Npv,Npv+2*Nref);
         Iref, sparse(Nref, N+Npv+2*Nref);
         sparse(Nref,N), Iref, sparse(Nref,Npv+2*Nref)];
     b = [P; Q; v2pv; 0; v2ref];
end

x = A\b;

if ~strcmp(version,'fatemi')
    theta = x(1:N);
    u     = x(N+1:2*N);
    v     = exp(u);
else
    thetap = x(1:N);
    v2     = x(N+1:2*N);
    theta  = theta_ref + thetap./v2;
    v      = sqrt(v2);
end

Qpv   = x(2*N+1:2*N+sum(pv_idx));
Pref  = x(2*N+sum(pv_idx)+1:2*N+sum(pv_idx) + sum(ref_idx));
Qref  = x(end-sum(ref_idx)+1:end);


If    = Yf*(v.*exp(1i*theta));
if ~strcmp(version,'fatemi')
    Pf    = v(nmap(mpc.branch(:,1))).*exp(1i*theta(nmap(mpc.branch(:,1)))).*conj(If);
else
    findx = sub2ind([N,N],nmap(mpc.branch(:,1)),nmap(mpc.branch(:,2)));
    Gtrf = Gtr(findx);
    Btrf = Btr(findx);
    E = sparse([1:M,1:M],[nmap(mpc.branch(:,1)); nmap(mpc.branch(:,2))],[ones(M,1);-ones(M,1)],M,N);
    Pf    = -Gtrf/2.*(E*thetap).^2 + Btrf.*(E*thetap);
    
    testp = zeros(N,1); testp(ref_idx) = Pref;
    testq = zeros(N,1); testq(ref_idx) = Qref;
    testq(pv_idx) = Qpv;
    test_ans = A(1:2*N,1:2*N)*[thetap;v2] - ([P;Q] + [testp; testq]);
end



%%
% findx = sub2ind([N,N],nmap(mpc.branch(:,1)),nmap(mpc.branch(:,2)));
% Gtrf = Gtr(findx);
% Btrf = Btr(findx);
% E = sparse([1:M,1:M],[nmap(mpc.branch(:,1)); nmap(mpc.branch(:,2))],[ones(M,1);-ones(M,1)],M,N);
% Sf = exp(abs(abs(E))*u).*(-Gtrf - Btrf.*E*theta + 1i*(-Gtrf.*E*theta + Btrf));
%%
mpopt = mpoption;
mpopt.out.all = 0;
mpcac = runpf(mpc,mpopt);
mpcac.If = Yf*(mpcac.bus(:,VM).*exp(1i*pi/180*mpcac.bus(:,VA)));
mpcac.Pf = mpcac.bus(nmap(mpcac.branch(:,1)),VM).*exp(1i*pi/180*mpcac.bus(nmap(mpcac.branch(:,1)),VA)).*conj(mpcac.If);
mpcdc = rundcpf(mpc,mpopt);

%% error
scag = struct('v',v - mpcac.bus(:,VM),...
    'theta', theta*180/pi - mpcac.bus(:,VA),...
    'Pf', real(Pf) - real(mpcac.Pf), 'Qf', imag(Pf) - imag(mpcac.Pf));
dc   = struct('v', 1-mpcac.bus(:,VM),...
    'theta', mpcdc.bus(:,VA) - mpcac.bus(:,VA),...
    'Pf', mpcdc.branch(:,PF)/mpc.baseMVA - real(mpcac.Pf));
%% case data
figure;
subplot(2,2,1);
plot(mpc.bus(:,PD));
hold on;
plot(mpc.bus(:,QD));
title('constant power load')

subplot(2,2,2)
plot(mpc.bus(:,GS));
hold on;
plot(mpc.bus(:,BS));
title('impedance load')

subplot(2,2,3)
plot(mpc.branch(:,BR_R));
hold on;
plot(mpc.branch(:,BR_X));
plot(mpc.branch(:,BR_B));
title('resistance,reactance, susceptance')

subplot(2,2,4)
plot(mpc.branch(:,TAP))
hold on;
plot(mpc.branch(:,SHIFT))
title('Tap and shift')

%%
figure;
plot(mpcac.bus(:,VM),'-o','linewidth',3);
hold on;
plot(v,'-*','linewidth',3)
plot(mpcdc.bus(:,VM),'-x','linewidth',3);
legend('ac pf', 'scaglione', 'dc pf')
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
plot(scag.v,'-o','linewidth',3);
ylabel('p.u. voltage error')
yyaxis right;
% plot(bsh/max(abs(bsh))*max(abs(scag.v)))
plot(bsh)
ylabel('shunt susceptance [p.u]')
xlabel('bus number')
ax = gca;
ax.FontSize = 16;

figure;
plot(scag.theta,'-o','linewidth',3);
hold on;
plot(dc.theta,'-x','linewidth',3);
legend('scaglione', 'dc pf')
xlabel('bus number')
ylabel('angle error [degree]')
ax = gca;
ax.FontSize = 16;

figure;
subplot(2,1,1)
plot(abs(mpcac.If),'-o','linewidth',3);
hold on;
plot(abs(If),'-x','linewidth',3);
legend('ac', 'scaglione')
xlabel('branch number')
ylabel('|I_f|')
ax = gca;
ax.FontSize = 16;
subplot(2,1,2)
plot(angle(mpcac.If)*180/pi,'-o','linewidth',3);
hold on;
plot(angle(If)*180/pi,'-x','linewidth',3);
legend('ac', 'scaglione')
xlabel('branch number')
ylabel('phase(I_f) [degree]')

figure;
subplot(2,1,1)
plot(abs(mpcac.Pf),'-o','linewidth',3);
hold on;
plot(abs(Pf),'-x','linewidth',3);
legend('ac', 'scaglione')
xlabel('branch number')
ylabel('|P_f|')
ax = gca;
ax.FontSize = 16;
subplot(2,1,2)
plot(angle(mpcac.Pf)*180/pi,'-o','linewidth',3);
hold on;
plot(angle(Pf)*180/pi,'-x','linewidth',3);
legend('ac', 'scaglione')
xlabel('branch number')
ylabel('phase(P_f) [degree]')

figure;
subplot(2,1,1)
plot(scag.Pf);
hold on;
plot(dc.Pf);
legend('scaglione', 'dc')
ylabel('real power error [p.u]')
ax = gca;
ax.FontSize = 16;
subplot(2,1,2)
plot(scag.Qf);
xlabel('branch number')
ylabel('reactive power error [p.u.]')
ax = gca;
ax.FontSize = 16;