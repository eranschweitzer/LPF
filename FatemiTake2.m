%%% Fatemi et al. take 2
clear variables; close all;
define_constants;
casename = 'case118';
mpc = loadcase(casename);
mpc0 = mpc;
lambda1 = 0.95;
lambda2 = 0.95;
%% remove out of service elements
mpc.bus    = mpc.bus(mpc.bus(:,BUS_TYPE) ~= 4,:);
mpc.branch = mpc.branch(mpc.branch(:,BR_STATUS) > 0, :);
mpc.gen    = mpc.gen(mpc.gen(:,GEN_STATUS) > 0, :);

% mpc.branch(:,TAP) = 0;
mpc.branch(:,BR_B) = 0;
% mpc.branch(:,SHIFT) = 0;
mpc.bus(:,GS) = 0;
mpc.bus(:,BS) = 0;
%%
nb = size(mpc.bus,1);
ng = size(mpc.gen,1);
nl = size(mpc.branch, 1);
%% sort case into REF | PV | PQ order
pqidx  = mpc.bus(:,BUS_TYPE) == PQ;
pvidx  = mpc.bus(:,BUS_TYPE) == PV;
refidx = mpc.bus(:,BUS_TYPE) == REF;
mx = 1:sum(refidx|pvidx);
nx = sum(refidx|pvidx)+1:nb;
mpc.bus = [mpc.bus(refidx,:); mpc.bus(pvidx,:); mpc.bus(pqidx,:)];
mpc.bus(1,VA) = 0;
%% renumber
e2i = sparse(mpc.bus(:,BUS_I), 1, 1:nb);
mpc.bus(:,BUS_I)    = full(e2i(mpc.bus(:,BUS_I)));
mpc.branch(:,F_BUS) = full(e2i(mpc.branch(:,F_BUS)));
mpc.branch(:,T_BUS) = full(e2i(mpc.branch(:,T_BUS)));
mpc.gen(:,GEN_BUS)  = full(e2i(mpc.gen(:,GEN_BUS)));

%% Load and generation
g2n = sparse(mpc.gen(:,GEN_BUS),1:ng,1,nb,ng);

Ptilde = (g2n*mpc.gen(:,PG) - mpc.bus(:,PD)) / mpc.baseMVA;
Qtilde = (g2n*mpc.gen(:,QG) - mpc.bus(:,QD)) / mpc.baseMVA;
P      = Ptilde(2:end);
Q      = Qtilde(2:end);
%% Y Matrix
[Y, Yf, Yt] = makeYbus(mpc);
G = lambda1*real(Y);
B = lambda2*imag(Y);
%% form H matrices
Bnninv = B(nx,nx)\speye(length(nx));
H11 = -(G(mx,nx)*Bnninv*G(nx,mx) + B(mx,mx)); %(B1)
H12 = -(G(mx,nx)*Bnninv*G(nx,nx) + B(mx,nx)); %(B2)
H21 = -(G(nx,nx)*Bnninv*G(nx,mx) + B(nx,mx)); %(B3)
H22 = -(G(nx,nx)*Bnninv*G(nx,nx) + B(nx,nx)); %(B4)

Htilde = [H11 H12; H21 H22]; %(39)
H = Htilde(2:end,2:end);

%% form L matrices
Lvm = G(mx,mx) - G(mx,nx)*Bnninv*B(nx,mx); %(B7)
LQm = -G(mx,nx)*Bnninv;                    %(B8)
Lvn = G(nx,mx) - G(nx,nx)*Bnninv*B(nx,mx); %(B9)
LQn = -G(nx,nx)*Bnninv;                    %(B10)

L   = [Lvm LQm; Lvn LQn];
%% make PvQ vector
vg = zeros(nb,1);
vg(mpc.gen(:,GEN_BUS)) = mpc.gen(:,VG);
v2m = vg(mx).^2;
Qn  = Qtilde(nx);
PvQtilde = L*[v2m;Qn]; %(40)
PvQ = PvQtilde(2:end);

%% solve
dp      = H\(P - PvQ); %(53)
dptilde = [mpc.bus(1,VA)*pi/180; dp];
v2n     = -Bnninv*(Qtilde(nx) + G(nx,mx)*dptilde(mx) + G(nx,nx)*dptilde(nx) + B(nx,mx)*v2m); %(35)
v       = sqrt([v2m;v2n]);
delta   = dptilde./[v2m;v2n];


%%
mpopt = mpoption;
mpopt.out.all = 0;
mpcac = runpf(mpc,mpopt);
% mpcac.If = Yf*(mpcac.bus(:,VM).*exp(1i*pi/180*mpcac.bus(:,VA)));
% mpcac.Pf = mpcac.bus(nmap(mpcac.branch(:,1)),VM).*exp(1i*pi/180*mpcac.bus(nmap(mpcac.branch(:,1)),VA)).*conj(mpcac.If);
mpcdc = rundcpf(mpc,mpopt);
%%
figure;
plot(mpcac.bus(:,VA),'-o','linewidth',3);
hold on;
plot(delta*180/pi,'-*','linewidth',3)
plot(mpcdc.bus(:,VA),'-x','linewidth',3);
legend('ac pf', 'Fatemi', 'dc pf')
xlabel('bus number')
ylabel('bus angle [degree]')
ax = gca;
ax.FontSize = 16;