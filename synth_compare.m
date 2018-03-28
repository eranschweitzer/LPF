function [Thist, fitness, Tsamp, C, mpc] = synth_compare(synmpc, casename, varargin)
%%% filename: name of synthetic case to be analyzed
%%% casename: name of Matpower case to compare to (optional)
%%% Name, Value:
%%%     - plot, (true/false): if true plots are generated. Default False

make_plots = varargin_parse(varargin,'plot', false);
plt_title  = varargin_parse(varargin,'title', '');
%% load case
define_constants;
mpopt = mpoption('opf.dc.solver', 'MIPS', 'opf.ac.solver', 'IPOPT', 'mips.step_control', 1);
% mpopt.opf.init_from_mpc = 1;
mpopt.out.all = 0;
mpopt.verbose = 0;
s = struct();
for prop = {'VMAX', 'RATE_A', 'PMAX', 'QMIN', 'QMAX'}
    %%% VMIN and PMIN are automatically set to with a bound of 0
    s.(prop{1}).type = 'unbnd';
end
for prop = {'ANGMIN', 'ANGMAX'}
    s.(prop{1}).type = 'none';
end
if ischar(synmpc)
    synmpc     = loadcase(synmpc);
end
synmpc.bus(:,VMAX) = 1.1;
synmpc.bus(:,VMIN) = 0.9;

synmpc.softlims = s;
synmpc   = toggle_softlims(synmpc,'on');
for k = 1:3
    mpcsolve = runopf(synmpc,mpopt);
    if mpcsolve.success
        break
    end
end
if ~mpcsolve.success
    Thist = NaN; fitness = NaN; Tsamp = NaN; C = NaN;
    warning('Synthetic case failed to converge')
    return
end
%%
N = size(synmpc.bus,1);
M = size(synmpc.branch,1);
nmap  = sparse(synmpc.bus(:,BUS_I),1,1:N);
E  = sparse([1:M,1:M].',[full(nmap(synmpc.branch(:,F_BUS)));full(nmap(synmpc.branch(:,T_BUS)))],[ones(M,1);-1*ones(M,1)],M,N);
F  = sparse(1:M,full(nmap(synmpc.branch(:,F_BUS))),1,M,N);
T  = sparse(1:M,full(nmap(synmpc.branch(:,T_BUS))),1,M,N);
Sb = BranchParts(synmpc);
vars.v    = synmpc.bus(:,VM);
vars.u0   = zeros(N,1);
vars.uhat = log(vars.v);
vars.theta= synmpc.bus(:,VA)*pi/180;
vars.phi  = 0.5*(E*vars.theta).^2;

P  = calcPflow(0,0,vars,F,T,E,Sb);
Q  = calcQflow(0,0,vars,F,T,E,Sb);
S  = struct('f', abs(P.f + 1i*Q.f), 't', abs(P.t + 1i*Q.t));
%% evaluate criteria between initial and AC solved case
Cv  = eval_criteria(synmpc.bus(:,VM),mpcsolve.bus(:,VM));
Ct  = eval_criteria(synmpc.bus(:,VA),mpcsolve.bus(:,VA));
Cpf = eval_criteria(P.f*synmpc.baseMVA,mpcsolve.branch(:,PF));
Cpt = eval_criteria(P.t*synmpc.baseMVA,mpcsolve.branch(:,PT));
Cqf = eval_criteria(Q.f*synmpc.baseMVA,mpcsolve.branch(:,QF));
Cqt = eval_criteria(Q.t*synmpc.baseMVA,mpcsolve.branch(:,QT));
Csf = eval_criteria(S.f*synmpc.baseMVA,abs(mpcsolve.branch(:,PF) + 1i*mpcsolve.branch(:,QF)));
Cst = eval_criteria(S.t*synmpc.baseMVA,abs(mpcsolve.branch(:,PT) + 1i*mpcsolve.branch(:,QT)));
Cdelta = eval_criteria(E*synmpc.bus(:,VA),E*mpcsolve.bus(:,VA));
C   = struct('v', Cv, 't', Ct, 'pf', Cpf, 'pt', Cpt, 'qf', Cqf, 'qt', Cqt,...
             'sf', Csf, 'st', Cst, 'delta', Cdelta);
%%
Tsamp = struct('v',synmpc.bus(:,VM),'vtrue',mpcsolve.bus(:,VM),...
       'delta',full(E*synmpc.bus(:,VA)),'deltatrue',full(E*mpcsolve.bus(:,VA)),...
       'pf', P.f*synmpc.baseMVA, 'pftrue', mpcsolve.branch(:,PF),...
       'qf', Q.f*synmpc.baseMVA, 'qftrue', mpcsolve.branch(:,QF),...
       'sf', S.f*synmpc.baseMVA, 'sftrue', abs(mpcsolve.branch(:,PF) + 1i*mpcsolve.branch(:,QF)));
Tsamp = mystruct2table(Tsamp);
%% comparison to a matpower case case
if ischar(casename)
    mpc = loadcase(casename);
    mpc.bus(:,VMAX) = 1.1;
    mpc.bus(:,VMIN) = 0.9;
    mpc.softlims = s;
    sgc = size(mpc.gencost,1);
    mpc.gencost(:,1) = 2;
    mpc.gencos = [mpc.gencost(:,1:3) 2*ones(sgc,1) 10*ones(sgc,1) zeros(sgc,1)];
    mpc = toggle_softlims(mpc, 'on');
    mpc = runopf(mpc,mpopt);
    if ~mpc.success
        error('Reference case failed to solve')
    end
else
    mpc = casename;
end
N2 = size(mpc.bus,1);
M2 = size(mpc.branch,1);
nmap2  = sparse(mpc.bus(:,BUS_I),1,1:N2);
E2  = sparse([1:M2,1:M2].',[full(nmap2(mpc.branch(:,F_BUS)));full(nmap2(mpc.branch(:,T_BUS)))],[ones(M2,1);-1*ones(M2,1)],M2,N2);

% T118 = struct('delta',full(E2*mpc.bus(:,VA)),'pf',mpc.branch(:,PF),...
%               'qf',mpc.branch(:,QF),'v',mpc.bus(:,VM));
%% Initialize histogram structure
Thist = struct();
Thist.q = (0.01:0.01:1).';
fitness = struct('kl', struct('delta', 0, 'pf', 0, 'qf', 0, 'sf', 0),...
                 'bc', struct('delta', 0, 'pf', 0, 'qf', 0, 'sf', 0));
%% Angle Difference
deltamax  = max([abs(Tsamp.deltatrue);abs(full(E2*mpc.bus(:,VA)))]);
Thist.deltaedges   = (-deltamax:0.5:deltamax+0.5);
Thist.deltacenters = Thist.deltaedges(1:end-1) + 0.5*diff(Thist.deltaedges);
[Thist.deltampc,~]   = histcounts(full(E2*mpc.bus(:,VA)),...
                                    Thist.deltaedges,'Normalization','pdf');
[Thist.deltasynth,~] = histcounts(Tsamp.deltatrue,...
                                    Thist.deltaedges,'Normalization','pdf');
Thist.deltampcq    = quantile(abs(full(E2*mpc.bus(:,VA))),Thist.q);
Thist.deltasynthq  = quantile(abs(Tsamp.deltatrue),Thist.q);

fitness.kl.delta = kld(Thist.deltasynth, Thist.deltampc, Thist.deltaedges, 'two_sample', true);
fitness.bc.delta = bc(Thist.deltasynth.*diff(Thist.deltaedges), Thist.deltampc.*diff(Thist.deltaedges));
%% Real Power (from end)
pfmax         = max([abs(Tsamp.pftrue);abs(mpc.branch(:,PF))]);
Thist.pfedges = (-pfmax:5:pfmax+5);
Thist.pfcenters = Thist.pfedges(1:end-1) + 0.5*diff(Thist.pfedges);
[Thist.pfmpc,~]   = histcounts(mpc.branch(:,PF),...
                               Thist.pfedges,'Normalization','pdf');
[Thist.pfsynth,~] = histcounts(Tsamp.pftrue,...
                               Thist.pfedges,'Normalization','pdf');
Thist.pfmpcq    = quantile(abs(mpc.branch(:,PF)),Thist.q);
Thist.pfsynthq  = quantile(abs(Tsamp.pftrue),Thist.q);

fitness.kl.pf = kld(Thist.pfsynth, Thist.pfmpc, Thist.pfedges, 'two_sample', true);
fitness.bc.pf = bc(Thist.pfsynth.*diff(Thist.pfedges), Thist.pfmpc.*diff(Thist.pfedges));
%% Reactive Power (from end)
qfmax         = max([abs(Tsamp.qftrue);abs(mpc.branch(:,QF))]);
Thist.qfedges = (-qfmax:5:qfmax+5);
Thist.qfcenters = Thist.qfedges(1:end-1) + 0.5*diff(Thist.qfedges);
[Thist.qfmpc,~]   = histcounts(mpc.branch(:,QF),...
                               Thist.qfedges,'Normalization','pdf');
[Thist.qfsynth,~] = histcounts(Tsamp.qftrue,...
                               Thist.qfedges,'Normalization','pdf');
Thist.qfmpcq    = quantile(abs(mpc.branch(:,QF)),Thist.q);
Thist.qfsynthq  = quantile(abs(Tsamp.qftrue),Thist.q);

fitness.kl.qf = kld(Thist.qfsynth, Thist.qfmpc, Thist.qfedges, 'two_sample', true);
fitness.bc.qf = bc(Thist.qfsynth.*diff(Thist.qfedges), Thist.qfmpc.*diff(Thist.qfedges));
%% Apparent Power (from-end)
Smpc   = mpc.branch(:,PF) + 1i*mpc.branch(:,QF);
Ssynth = mpcsolve.branch(:,PF) + 1i*mpcsolve.branch(:,QF);
sfmax         = max([abs(Ssynth);abs(Smpc)]);
Thist.sfedges = (-sfmax:5:sfmax+5);
Thist.sfcenters = Thist.sfedges(1:end-1) + 0.5*diff(Thist.sfedges);
[Thist.sfmpc,~]   = histcounts(abs(Smpc),...
                               Thist.sfedges,'Normalization','pdf');
[Thist.sfsynth,~] = histcounts(abs(Ssynth),...
                               Thist.sfedges,'Normalization','pdf');
Thist.sfmpcq    = quantile(abs(Smpc),Thist.q);
Thist.sfsynthq  = quantile(abs(Ssynth),Thist.q);

fitness.kl.sf = kld(Thist.sfsynth, Thist.sfmpc, Thist.sfedges, 'two_sample', true);
fitness.bc.sf = bc(Thist.sfsynth.*diff(Thist.sfedges), Thist.sfmpc.*diff(Thist.sfedges));

%% SIL
SILmpc   = SIL_loading(mpc, 'xfmr_chk', 'lax');
SILsynth = SIL_loading(mpcsolve, 'xfmr_chk', 'lax');
silmax = max([SILmpc; SILsynth]);
Thist.siledges   = (0:0.1:silmax+0.1);
Thist.silcenters = Thist.siledges(1:end-1) + 0.5*diff(Thist.siledges);
[Thist.silmpc,~]   = histcounts(SILmpc,...
                               Thist.siledges,'Normalization','pdf');
[Thist.silsynth,~] = histcounts(SILsynth,...
                               Thist.siledges,'Normalization','pdf');
Thist.silmpcq    = quantile(SILmpc,Thist.q);
Thist.silsynthq  = quantile(SILsynth,Thist.q);

fitness.kl.sil = kld(Thist.silsynth, Thist.silmpc, Thist.siledges, 'two_sample', true);
fitness.bc.sil = bc(Thist.silsynth.*diff(Thist.siledges), Thist.silmpc.*diff(Thist.siledges));

mu = mean(SILsynth);
fitness.sil_exp.synth.mu = mu;
fitness.sil_exp.synth.kl = kld(Thist.silsynth, 1/mu*exp(-Thist.silcenters/mu), Thist.siledges);
fitness.sil_exp.synth.bc = bc(Thist.silsynth.*diff(Thist.siledges), 1/mu*exp(-Thist.silcenters/mu).*diff(Thist.siledges));
mu = mean(SILmpc);
fitness.sil_exp.real.mu  = mu;  
fitness.sil_exp.real.kl = kld(Thist.silmpc, 1/mu*exp(-Thist.silcenters/mu), Thist.siledges);
fitness.sil_exp.real.bc = bc(Thist.silmpc.*diff(Thist.siledges), 1/mu*exp(-Thist.silcenters/mu).*diff(Thist.siledges));
%% Figures
if ~make_plots
    return
end
synth_compare_plots(Thist,'title', plt_title)
% % Angle
% figure;
% subplot(1,2,1)
% plot(Thist.deltacenters,Thist.deltampc,'-o');
% hold on;
% plot(Thist.deltacenters,Thist.deltasynth,'-*');
% set(gca,'Yscale','log')
% xlabel('Angel Difference [degrees]')
% ylabel('Density')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% 
% subplot(1,2,2)
% plot(Thist.q,Thist.deltampcq);
% hold on;
% plot(Thist.q,Thist.deltasynthq);
% set(gca,'Yscale','log')
% ylabel('Angel Difference [degrees]')
% xlabel('quantile')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% ax = gca;
% ax.YLim(1) = 1e-5;
% 
% % Real Powerflow (from end)
% figure;
% subplot(1,2,1)
% plot(Thist.pfcenters,Thist.pfmpc,'-o');
% hold on;
% plot(Thist.pfcenters,Thist.pfsynth,'-*');
% set(gca,'Yscale','log')
% xlabel('Real Power (from-end) [MW]')
% ylabel('Density')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% 
% subplot(1,2,2)
% plot(Thist.q,Thist.pfmpcq);
% hold on;
% plot(Thist.q,Thist.pfsynthq);
% set(gca,'Yscale','log')
% ylabel('Real Power (from-end) [MW]')
% xlabel('quantile')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% ax = gca;
% ax.YLim(1) = 1e-5;
% 
% % Reactive Powerflow (from end)
% figure;
% subplot(1,2,1)
% plot(Thist.qfcenters,Thist.qfmpc,'-o');
% hold on;
% plot(Thist.qfcenters,Thist.qfsynth,'-*');
% set(gca,'Yscale','log')
% xlabel('Reactive Power (from-end) [MVAr]')
% ylabel('Density')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% 
% subplot(1,2,2)
% plot(Thist.q,Thist.qfmpcq);
% hold on;
% plot(Thist.q,Thist.qfsynthq);
% set(gca,'Yscale','log')
% ylabel('Reactive Power (from-end) [MVAr]')
% xlabel('quantile')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% ax = gca;
% ax.YLim(1) = 1e-5;
% 
% % Apparent Powerflow (from end)
% figure;
% subplot(1,2,1)
% plot(Thist.sfcenters,Thist.sfmpc,'-o');
% hold on;
% plot(Thist.sfcenters,Thist.sfsynth,'-*');
% set(gca,'Yscale','log')
% xlabel('Apparent Power (from-end) [MVA]')
% ylabel('Density')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% 
% subplot(1,2,2)
% plot(Thist.q,Thist.sfmpcq);
% hold on;
% plot(Thist.q,Thist.sfsynthq);
% set(gca,'Yscale','log')
% ylabel('Apparent Power (from-end) [MVA]')
% xlabel('quantile')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% ax = gca;
% ax.YLim(1) = 1e-5;