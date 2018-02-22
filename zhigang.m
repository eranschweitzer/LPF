function vars = zhigang(mpc)
% Direct implementation of Zhigang's paper
%%%%%% THERE IS AN ERROR SOMEWHERE %%%%%%%%%%%%%%%
define_constants;
if ischar(mpc)
    mpc = loadcase(mpc);
end

N = struct('t',size(mpc.bus,1));
M = size(mpc.branch,1);
G = size(mpc.gen,1);
baseMVA = mpc.baseMVA;
%% maps and operating matrices
nmap  = sparse(mpc.bus(:,BUS_I),1,1:N.t);
E  = sparse([1:M,1:M].',[full(nmap(mpc.branch(:,F_BUS)));full(nmap(mpc.branch(:,T_BUS)))],[ones(M,1);-1*ones(M,1)],M,N.t);
F  = sparse(1:M,full(nmap(mpc.branch(:,F_BUS))),1,M,N.t);
T  = sparse(1:M,full(nmap(mpc.branch(:,T_BUS))),1,M,N.t);
gbus  = full(nmap(mpc.gen(:,GEN_BUS))); %generator buses
gstat = mpc.gen(:,GEN_STATUS) > 0;  % mask of dispatched generators
gmap  = sparse(gbus,1:G,gstat,N.t,G); %maps generators onto buses
bus_with_ongen = sum(gmap,2) > 0;

%% bus types
bidx = bustype_map(mpc,bus_with_ongen,gstat,nmap);
N.pq = sum(bidx.pq); N.pv = sum(bidx.pv); N.ref= sum(bidx.ref);
%% u0
vg = ones(N.t,1); vg(gbus(gstat)) = mpc.gen(gstat,VG); 
u0 = u0init(vg,bidx,'udefault',0);
theta_ref = mpc.bus(bidx.ref,VA)*pi/180;

%% load and generation
Sp = powervectors(mpc,gmap);

%% branch parts
Sb = BranchParts(mpc);

%% real part
% equation (19)
P = (Sp.Pg - Sp.Pd) - Sb.gsh - F'*(-(Sb.g.*log(Sb.tau) - Sb.b.*Sb.tshift)./Sb.tau) ...
    -T'*(Sb.g.*log(Sb.tau) - Sb.b.*Sb.tshift);

% equation (20)
gshp = (Sp.Pg - Sp.Pd) + Sb.gsh;

% equation (21)
gp   = F'*diags(Sb.g./Sb.tau)*E - T'*diags(Sb.g)*E;
bp   = F'*diags(Sb.b./Sb.tau)*E - T'*diags(Sb.b)*E;

%% reactive part
% equation (28)
Q = (Sp.Qd - Sp.Qg) - Sb.bsh - F'*( (Sb.bc/2 - (Sb.b.*log(Sb.tau) + Sb.g.*Sb.tshift))./Sb.tau.^2)...
    -T'*(Sb.bc/2 + Sb.b.*log(Sb.tau) + Sb.g.*Sb.tshift);

% equation (29)
gpp = F'*diags(Sb.g./Sb.tau.^2)*E - T'*diags(Sb.g)*E;
bpp = F'*diags(Sb.b./Sb.tau.^2)*E - T'*diags(Sb.b)*E;

%%
Art = -bp;
Aru = diags(gshp) + gp;
Amt = gpp;
Amu = diags(2*(Sp.Qd - Sp.Qg)) + bpp;

matI = struct(); 
matI.pv  = sparse(1:N.pv,find(bidx.pv),1, N.pv, N.t);
matI.ref  = sparse(1:N.ref,find(bidx.ref),1, N.ref, N.t);

matU = struct();
matU.pv  = sparse(find(bidx.pv),1:N.pv,1-2*u0(bidx.pv), N.t, N.pv);
matU.ref = sparse(find(bidx.ref),1:N.ref,1-2*u0(bidx.ref), N.t, N.ref);
matU.refp = sparse(find(bidx.ref),1:N.ref,-1 + u0(bidx.ref), N.t, N.ref);

A = [ Art              , Aru              , sparse(N.t,N.pv), matU.refp          , sparse(N.t, N.ref);
      Amt              , Amu              , matU.pv         , sparse(N.t, N.ref) , matU.ref;
      sparse(N.pv,N.t) , matI.pv          , sparse(N.pv, N.pv + 2*N.ref);
      matI.ref         , sparse(N.ref,N.t), sparse(N.ref, N.pv + 2*N.ref);
      sparse(N.ref,N.t), matI.ref         , sparse(N.ref, N.pv + 2*N.ref)];
b = [P;Q; u0(bidx.pv); theta_ref; u0(bidx.ref)];

x = A\b;
vars   = result_parse(x,u0,N);
vars.v = exp(vars.uhat);
vars.u0(:) = 0;
flows = zhigangflow(vars,F,T,E,Sb);
vars.sf = flows.f;
vars.st = flows.t;