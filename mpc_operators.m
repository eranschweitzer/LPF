function S = mpc_operators(mpc)
%%% useful operator matrices and maps when working with the mpc structure

define_constants;
%% sizes
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
%% form structure
S = struct('N',N, 'M', M, 'G', G, 'baseMVA', baseMVA, 'nmap', nmap,...
    'E', E, 'F', F, 'T', T, 'gbus', gbus, 'gstat', gstat, 'gmap', gmap, ...
    'bus_with_ongen', bus_with_ongen, 'bidx', bidx);