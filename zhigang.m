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

%% real part (not that re(y^sh) = 0 i.e. only branch shunt susceptance)
% equation (19)
Pp = (Sp.Pg - Sp.Pd) - Sb.gsh - F'*(-(Sb.g.*log(Sb.tau) - Sb.b.*Sb.tshift)./Sb.tau) ...
    -T'*(Sb.g.*log(Sb.tau) - Sb.b.*Sb.tshift);

% equation (20)
gshp = (Sp.Pg - Sp.Pd) + Sb.gsh;

% equation (23)
% gpold   = F'*diags(Sb.g./Sb.tau)*E - T'*diags(Sb.g)*E;
gp   = -(F'*diags(Sb.g./Sb.tau)*T + T'*diags(Sb.g)*F);
gpt  = diags(gshp) + (F'*diags(Sb.g./Sb.tau)*F + T'*diags(Sb.g)*T);
% gpt  = gpt + gp;
% equation (24)
bp   = F'*diags(Sb.b./Sb.tau)*E - T'*diags(Sb.b)*E;

%% reactive part
% equation (28)
Qpp = -( (Sp.Qg - Sp.Qd) - Sb.bsh - F'*( (Sb.bc/2 - (Sb.b.*log(Sb.tau) + Sb.g.*Sb.tshift))./Sb.tau.^2)...
    -T'*(Sb.bc/2 + Sb.b.*log(Sb.tau) + Sb.g.*Sb.tshift) );

% equation (29)
gpp  = F'*diags(Sb.g./Sb.tau.^2)*E - T'*diags(Sb.g)*E;
% bpp = F'*diags(Sb.b./Sb.tau.^2)*E - T'*diags(Sb.b)*E;
bpp  = -(F'*diags(Sb.b./Sb.tau.^2)*T + T'*diags(Sb.b)*F);
bppt = -2*diags(Sp.Qg - Sp.Qd) + (F'*diags(Sb.b./Sb.tau.^2)*F + T'*diags(Sb.b)*T);
% bppt = bppt + bpp;
%% equation (53)
A = [bp(bidx.pq, bidx.ref)  -gpt(bidx.pq,bidx.pv)  -gpt(bidx.pq,bidx.ref);
     bp(bidx.pv, bidx.ref)  -gpt(bidx.pv,bidx.pv)  -gpt(bidx.pv,bidx.ref);
     gpp(bidx.pq, bidx.ref)  bppt(bidx.pq,bidx.pv)  bppt(bidx.pq,bidx.ref)];

vtilde = [Pp(bidx.pq); Pp(bidx.pv); Qpp(bidx.pq)] + ...
          A*[theta_ref; u0(bidx.pv); u0(bidx.ref)];
 
%% equation(52)
A = [bp(bidx.pq,bidx.pq)  bp(bidx.pq,bidx.pv) -gpt(bidx.pq,bidx.pq);
     bp(bidx.pv,bidx.pq)  bp(bidx.pv,bidx.pv) -gpt(bidx.pv,bidx.pq);
     gpp(bidx.pq,bidx.pq) gpp(bidx.pq,bidx.pv) bppt(bidx.pq,bidx.pq)];

x = -A\vtilde;

%% variables
vars = struct();
vars.theta = zeros(N.t,1);
vars.u     = zeros(N.t,1);

% ref
vars.theta(bidx.ref) = theta_ref;
vars.u(bidx.ref)     = u0(bidx.ref);

% pv
vars.theta(bidx.pv) = x(N.pq+1:N.pq+N.pv);
vars.u(bidx.pv)     = u0(bidx.pv);

%pq
vars.theta(bidx.pq) = x(1:N.pq);
vars.u(bidx.pq)     = x(N.pq+N.pv+1:end);
% %%
% Art = -bp;
% Aru = diags(gshp) + gp;
% Amt = gpp;
% Amu = diags(2*(Sp.Qd - Sp.Qg)) + bpp;
% 
% matI = struct(); 
% matI.pv  = sparse(1:N.pv,find(bidx.pv),1, N.pv, N.t);
% matI.ref  = sparse(1:N.ref,find(bidx.ref),1, N.ref, N.t);
% 
% matU = struct();
% matU.pv  = sparse(find(bidx.pv),1:N.pv,1-2*u0(bidx.pv), N.t, N.pv);
% matU.ref = sparse(find(bidx.ref),1:N.ref,1-2*u0(bidx.ref), N.t, N.ref);
% matU.refp = sparse(find(bidx.ref),1:N.ref,-1 + u0(bidx.ref), N.t, N.ref);
% 
% A = [ Art              , Aru              , sparse(N.t,N.pv), matU.refp          , sparse(N.t, N.ref);
%       Amt              , Amu              , matU.pv         , sparse(N.t, N.ref) , matU.ref;
%       sparse(N.pv,N.t) , matI.pv          , sparse(N.pv, N.pv + 2*N.ref);
%       matI.ref         , sparse(N.ref,N.t), sparse(N.ref, N.pv + 2*N.ref);
%       sparse(N.ref,N.t), matI.ref         , sparse(N.ref, N.pv + 2*N.ref)];
% b = [P;Q; u0(bidx.pv); theta_ref; u0(bidx.ref)];
% 
% x = A\b;
% vars   = result_parse(x,u0,N);
vars.v = exp(vars.u);
% vars.u0(:) = 0;
flows = zhigangflow(vars,F,T,E,Sb);
vars.sf = flows.f;
vars.st = flows.t;