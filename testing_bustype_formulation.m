clear variables; %close all;
%% inputs
casename = 'case3375wp';%'case1888rte';%'case_ACTIVSg2000';%'case24_ieee_rts';%'case3375wp';%'case1354pegase';%'case6470rte';%'case13659pegase';%'case_ACTIVSg2000';%%;%'case6470rte'%'case_ACTIVSg10k';
%% initializations
% define_constatns;
mpc = loadcase(casename);  
opmats = mpc_operators(mpc);
v2struct(opmats);
%% branch parts
Sb = BranchParts(mpc);
%% power vectors
Sp = powervectors(mpc,gmap);
%% vtrue
vtrue = vtrue_struct(mpc,opmats, Sb, Sp);
%% u0
vg = ones(N.t,1); 
vg(gbus(gstat)) = mpc.gen(gstat,6); % 6=VG in gen matrix
% unearest = nearestu0(F,T,Sb,N,bidx);
% u0 = u0init(vg,bidx,'udefault','nearest','unearest',unearest);
u0 = u0init(vg,bidx,'udefault',0);
% u0 = log(vtrue.v);
tref = mpc.bus(bidx.ref,9)*pi/180; %9=VA in bus matrix
%% reweight Y
% yw = reweightY(opmats,Sb,Sp,u0,'pqopt',1);

%% A matrix
% b = struct();
% A = struct('pq', struct(), 'pv', struct());
% b.pq = conj(yw.S(bidx.pq)) - (yw.ysh(bidx.pq) + yw.yff(bidx.pq) + yw.ytt(bidx.pq) + ...
%             cnx.F.pq'*(yw.yft + 1i*tref*(yw.yft.*lidx.T.ref)) + ...
%             cnx.T.pq'*(yw.ytf + 1i*tref*(yw.ytf.*lidx.F.ref)));
% A.pq.tpq = 1i*( (-cnx.F.pq'*diags(yw.yft) + cnx.T.pq'*diags(yw.ytf.*lidx.F.pq))*cnx.F.pq +...
%                 (+cnx.F.pq'*diags(yw.yft.*lidx.T.pq) - cnx.T.pq'*diags(yw.ytf))*cnx.T.pq);
% A.pq.tpv = 1i*(cnx.F.pq'*diags(yw.yft.*lidx.T.pv)*cnx.T.pv + cnx.T.pq'*diags(yw.ytf.*lidx.F.pv)*cnx.F.pv);
% switch yw.pqopt
%     case 1
%         A.pq.uhat = diags(conj(yw.S(bidx.pq)) + yw.ysh(bidx.pq) + yw.yff(bidx.pq) + yw.ytt(bidx.pq)) + ...
%                     cnx.F.pq'*diags(yw.yft.*lidx.T.pq)*cnx.T.pq + ...
%                     cnx.T.pq'*diags(yw.ytf.*lidx.F.pq)*cnx.F.pq;
%     case 2
%         A.pq.uhat = A.pq.tpq/1i;
%     otherwise
%         error('yw.pqopt must be either 1 or 2')
% end
% 
% b.pv = conj(yw.S(bidx.pv)) - (yw.ysh(bidx.pv) + yw.yff(bidx.pv) + yw.ytt(bidx.pv) + ...
%             cnx.F.pv'*(yw.yft + 1i*tref*(yw.yft.*lidx.T.ref)) + ...
%             cnx.T.pv'*(yw.ytf + 1i*tref*(yw.ytf.*lidx.F.ref)));
% A.pv.tpq = 1i*(cnx.F.pv'*diags(yw.yft.*lidx.T.pq)*cnx.T.pq + cnx.T.pv'*diags(yw.ytf.*lidx.F.pq)*cnx.F.pq);
% A.pv.uhat= A.pv.tpq/1i;
% A.pv.tpv = 1i*( (-cnx.F.pv'*diags(yw.yft) + cnx.T.pv'*diags(yw.ytf.*lidx.F.pv))*cnx.F.pv +...
%                 (+cnx.F.pv'*diags(yw.yft.*lidx.T.pv) - cnx.T.pv'*diags(yw.ytf))*cnx.T.pv);
% A.pv.qg  = 1i*speye(N.pv);

%% solve
% Amat = [real(A.pq.tpq), real(A.pq.uhat), real(A.pq.tpv), sparse(N.pq,N.pv);
%      imag(A.pq.tpq), imag(A.pq.uhat), imag(A.pq.tpv), sparse(N.pq,N.pv);
%      real(A.pv.tpq), real(A.pv.uhat), real(A.pv.tpv), real(A.pv.qg);
%      imag(A.pv.tpq), imag(A.pv.uhat), imag(A.pv.tpv), imag(A.pv.qg)];
% b = [real(b.pq); imag(b.pq); real(b.pv); imag(b.pv)] + ...
%     [real(A.pq.uhat); imag(A.pq.uhat); real(A.pv.uhat); imag(A.pv.uhat)]*u0(bidx.pq);
for k = 1

if k > 1
    phi = ((F-T)*vars.theta).^2/2;
%     phi(phi > (30*pi/180)^2/2) = 0;
else
    phi = 0*((F-T)*vtrue.t).^2/2;% - 1i*((F-T)*vtrue.t).^2/6;
end

if k == 3
    u0(bidx.pq) = vars.u(bidx.pq);
end

[Amat, bmat, Apmat] = blockMats(opmats,Sb,Sp,u0,tref);
phitmp = [Apmat{1}+1i*Apmat{2}; Apmat{3}+1i*Apmat{4}]*phi;
b = cell2mat(bmat) + cell2mat(Amat(:,2))*u0(bidx.pq) - ...
    [real(phitmp(1:N.pq)); imag(phitmp(1:N.pq)); real(phitmp(N.pq+1:end)); imag(phitmp(N.pq+1:end))];
% b = cell2mat(bmat) + cell2mat(Amat(:,2))*u0(bidx.pq) - cell2mat(Apmat)*phi;
A = cell2mat(Amat);
% x =(0.01*speye(size(Amat)) + Amat'*Amat)\(Amat'*b);
x = A\b;

vars = struct('u', u0, 'theta',zeros(N.t,1));
vars.theta(bidx.ref)= tref;
vars.theta(bidx.pq) = x(1:N.pq);
vars.theta(bidx.pv) = x(2*N.pq+1:2*N.pq+N.pv);
vars.u(bidx.pq)     = x(N.pq+1:2*N.pq);
vars.Qpv            = x(end-N.pv+1:end);

vars.v = exp(vars.u + 1i*vars.theta);

Ybus = myMakeYbus(F,T,Sb);
sref = vars.v(bidx.ref)*conj(Ybus(bidx.ref,:)*vars.v) + Sp.Pd(bidx.ref) + 1i*Sp.Qd(bidx.ref);
vars.Pref = real(sref);
vars.Qref = imag(sref);

Sg   = makeSg(vars,Sp,bidx,N);
residual = pfresidual(vars.v,Ybus,real(Sg) - Sp.Pd, imag(Sg) - Sp.Qd);
%% plots
vt_comp_plots(vars, vtrue, E, residual,bidx)
end