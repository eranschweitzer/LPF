function vars = SLPF(mpc,varargin)

%% optional inputs
itermax  = varargin_parse(varargin,'itermax',50);
GwAru    = varargin_parse(varargin,'GwAru','default');
GwAmu    = varargin_parse(varargin,'GWAmu','default');
bcmode   = varargin_parse(varargin,'bcmode','default');
bshmode  = varargin_parse(varargin,'bshmode','default');
udefault = varargin_parse(varargin,'udefault','uref');
uopt     = varargin_parse(varargin,'uopt','none');
% id       = varargin_parse(varargin,'id','0315010010110');
id       = varargin_parse(varargin,'id','022200000000');
%% load case
define_constants;
if ischar(mpc)
    mpc = loadcase(mpc);
elseif ~isstruct(mpc)
    error('mpc must be either a struct or a casename string')
end

N = struct('t',size(mpc.bus,1));
M = size(mpc.branch,1);
G = size(mpc.gen,1);
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
u0 = u0init(vg,bidx,'udefault',udefault);
theta_ref = mpc.bus(bidx.ref,VA)*pi/180;

%% load and generation
Sp = powervectors(mpc,gmap);

%% branch parts
Sb = BranchParts(mpc);
Gw = branchweights(Sb,'GwAru',GwAru,'GwAmu',GwAmu,'bcmode',bcmode,'bshmode',bshmode);
Sb.b    = Gw.b*Sb.b;
Sb.balt = Gw.b*Sb.balt;
Sb.bc   = Gw.bc*Sb.bc;
Sb.bsh  = Gw.bsh*Sb.bsh;
Sb.g    = Gw.g*Sb.g;
Gw = branchweights(Sb);

%% Main Loop

ids  = str2ids(id);
vars = pfsolve(ids,F,T,E,Sb,Sp,bidx,theta_ref,N,u0,...
               'itermax',itermax,'Gw',Gw,'uopt',uopt);
if strcmp(id,'022200000000')
  vars.p = zhigangflow(vars, F, T, E, Sb, 'pq', 'real');
  vars.q = zhigangflow(vars, F, T, E, Sb, 'pq', 'imag');
end