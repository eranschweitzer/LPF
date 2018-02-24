clear variables; close all;
%% load case
define_constants;
% casename = '/Users/eran/Dropbox/ASU/SINE/python/parameter_assignment/test118minloss';
casename = 'case1888rte';%'case118';%'case24_ieee_rts';%'case3375wp';%'case1354pegase';%'case6470rte';%'case13659pegase';%'case_ACTIVSg2000';%%;%'case6470rte'%'case_ACTIVSg10k';
mpc = loadcase(casename);
% id = '0000000000000';
% id = '011000000000';
% id = '022220111001';
id = '022200000000';
% id = '1000000000000';
% id  = '0315010010110';
% id  = '0005000100000';
% id  = '3005010000111';
% id  = '0006000000000';
itermax = 75;
GwAru = 'default';
GwAmu = 'default';
bcmode= 'default';
bshmode= 'default';
udefault = 'uref';%{0,'uref'};
uopt   = 'none';
lsqr   = false;
plots  = true;
id = '022020011000'; %apparent winner
% id = {'011021011010', '021021011010'};%, '021220011111'};
% id = '021220011111';

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
u0 = u0init(vg,bidx,'udefault',udefault);
theta_ref = mpc.bus(bidx.ref,VA)*pi/180;

%% matrices fixing PV and Ref quantities
% [matI,matU] = pv_ref_mats(bidx,N,u0);
%% load and generation
Sp = powervectors(mpc,gmap);
Sp.Qg(bidx.pv) = 0;
Sp.Qg(bidx.ref) = 0;
Sp.Pg(bidx.ref) = 0;
%% branch parts
Sb = BranchParts(mpc);
Gw = branchweights(Sb,'GwAru',GwAru,'GwAmu',GwAmu,'bcmode',bcmode,'bshmode',bshmode);
Sb.b    = Gw.b*Sb.b;
Sb.balt = Gw.b*Sb.balt;
Sb.bc   = Gw.bc*Sb.bc;
Sb.bsh  = Gw.bsh*Sb.bsh;
Sb.g    = Gw.g*Sb.g;
Gw = branchweights(Sb);
%% AC solution
mpopt   = mpoption; mpopt.out.all = 0;
mpcac   = runpf(mpc,mpopt);
vtrue   = struct();
vtrue.v = mpcac.bus(:,VM);        %true voltage magnitude solution
vtrue.t = mpcac.bus(:,VA)*pi/180; %true voltage angle solution
vtrue.phi = 0.5*(E*vtrue.t).^2;       %true phi (1/2)*(delta(theta))^2
vtrue.pf  = mpcac.branch(:,PF);     %true real power at from bus [MW]
vtrue.pt  = mpcac.branch(:,PT);     %true real power at to bus [MW]
vtrue.qf  = mpcac.branch(:,QF);     %true reactive power at from bus [MVAr]
vtrue.qt  = mpcac.branch(:,QT);     %true reactive power at to bus [MVAr]
%% Main Loop

ids  = str2ids(id);
if lsqr
    ids = makeop(ids,u0);
end
vars = pfsolve(ids,F,T,E,Sb,Sp,bidx,theta_ref,N,u0,...
               'itermax',itermax,'Gw',Gw,'uopt',uopt, 'lsqr', lsqr);
if ~vars.convg
    return
end
% Evaluation criteria for voltage and angle 
Cv   = eval_criteria(vars.v,vtrue.v);
Ct   = eval_criteria(vars.theta,vtrue.t);
if plots
    vt_comp_plots(vars,vtrue)
end
% Evalutaion criteria for flow calculation
if length(id) ~= 13
    flids = flowids();
%     Cpf = cell(length(flids.p),1); Cpt = cell(length(flids.p),1);
%     Cqf = cell(length(flids.q),1); Cqt = cell(length(flids.q),1);
%     for k = 1:length(flids.p)
%         flid = str2ids(flids.p{k});
%         flid.c = ids.c;
%         if flid.Q == 3
%         end
%         P = flowcalc(flid, 'real', vars, F, T, E, Sb);
%         Cpf{k} = eval_criteria(P.f, Pftrue/baseMVA);
%         Cpt{k} = eval_criteria(P.t, Pttrue/baseMVA);
%     end
%     for k = 1:length(flids.q)
%         flid = str2ids(flids.q{k});
%         flid.c = ids.c;
%         if flid.Q == 3
%         end
%         Q = flowcalc(flid, 'imag', vars, F, T, E, Sb);
%         Cqf{k} = eval_criteria(Q.f, Qftrue/baseMVA);
%         Cqt{k} = eval_criteria(Q.t, Qttrue/baseMVA);
%     end
    Cf = flow_performance(flids, ids, vars, vtrue, F, T, E, Sb, baseMVA);
    return
else
    Cpf = cell(6,2); Cpt = cell(6,2);
    Cqf = cell(6,2); Cqt = cell(6,2);
    for ID = 0:5
        for btype = 0:1
            P = calcPflow(ID,btype,vars,F,T,E,Sb,'Gw',Gw);
            Q = calcQflow(ID,btype,vars,F,T,E,Sb,'Gw',Gw);
            Cpf{ID+1,btype+1} = eval_criteria(baseMVA*P.f,vtrue.pf);
            Cpt{ID+1,btype+1} = eval_criteria(baseMVA*P.t,vtrue.pt);
            Cqf{ID+1,btype+1} = eval_criteria(baseMVA*Q.f,vtrue.qf);
            Cqt{ID+1,btype+1} = eval_criteria(baseMVA*Q.t,vtrue.qt);
        end
    end

    % create result table 
    colnames = {'casename', 'id', 'prop', 'criteria', 'value'};
    Tres = result2table(colnames,casename,ids,vars.convg,Cv,Ct,Cpf,Cpt,Cqf,Cqt);
end
%% testing vecotrs
vect = struct();
vect.max = strcmp(Tres.criteria,'max');
vect.avg = strcmp(Tres.criteria,'avg');
vect.cor = strcmp(Tres.criteria,'cor');
vect.del = strcmp(Tres.criteria,'del');
%%
Q = calcQflow(0,0,vars,F,T,E,Sb,'Gw',Gw);
[~,idx] = max(abs(Qftrue - Q.f*baseMVA));
fnode = find(E(idx,:) == 1);
tnode = find(E(idx,:) == -1);
shunt = struct('from',Sb.bsh(fnode),'to',Sb.bsh(tnode));
bchar = struct('from',F(:,fnode)'*Sb.bc,'to',T(:,tnode)'*Sb.bc);

Qfparts       = struct;
Qfparts.fu    = -exp(Gw.b*F*vars.u0).*Sb.tau.^(-2).*(Sb.b + Sb.bc/2).*(1 + F*vars.uhat);
Qfparts.tu    = +exp(Gw.b*T*vars.u0).*Sb.tau.^(-1).*(Sb.g.*sin(Sb.tshift) + Sb.b.*cos(Sb.tshift)).*(1 + T*vars.uhat - vars.phi);
Qfparts.theta = -exp(T*vars.u0).*Sb.tau.^(-1).*(Sb.g.*cos(Sb.tshift) - Sb.b.*sin(Sb.tshift)).*E*vars.theta;
Qfparts.p     = F*vars.v;
Qfparts.full  = Qfparts.p.*(Qfparts.fu + Qfparts.tu + Qfparts.theta);

Qrparts       = struct;
Qrparts.fu    = -exp(Gw.b*F*log(vtrue)).*Sb.tau.^(-2).*(Sb.b + Sb.bc/2);
Qrparts.tu    = +exp(Gw.b*T*log(vtrue)).*Sb.tau.^(-1).*(Sb.g.*sin(Sb.tshift) + Sb.b.*cos(Sb.tshift)).*(1 - phitrue);
Qrparts.theta = -exp(T*log(vtrue)).*Sb.tau.^(-1).*(Sb.g.*cos(Sb.tshift) - Sb.b.*sin(Sb.tshift)).*E*ttrue;
Qrparts.p     = F*vtrue;
Qrparts.full  = Qrparts.p.*(Qrparts.fu + Qrparts.tu + Qrparts.theta);

%% linpf start vs flat
actest = linpf2ac_convergence(mpc,vars,bidx);