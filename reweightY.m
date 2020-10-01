function out = reweightY(opmats,Sb,Sp,u0,varargin)

pqopt = varargin_parse(varargin,'pqopt',1);
gs = varargin_parse(varargin,'gs',0);
rsimp  = varargin_parse(varargin,'rsimp',0);
isimp  = varargin_parse(varargin,'isimp',0);
vtau   = varargin_parse(varargin,'vtau',true);
v2struct(opmats); %move fields to variables
%% branch parts
[g,b,bc,tau,tshift,gsh,bsh] = branchparts_unpack(Sb,gs>0);
if gs > 0
    g = zeros(size(g));
end
if rsimp == 1
    % simplification for real power is setting tau=1
    tau = ones(size(tau));
end
if isimp == 1
    % simplification for imaginary/reactive power means setting tshift = 0
    tshift = zeros(size(tshift));
end
%% define base values
ys = g + 1i*b;
ytt = ys + 0.5*1i*bc;
yff = (ytt)./tau.^2;
yft = -ys.*(exp(+1i*tshift)./tau);
ytf = -ys.*(exp(-1i*tshift)./tau);
ysh = gsh + 1i*bsh;
S   = (Sp.Pg - Sp.Pd) + 1i*(Sp.Qg - Sp.Qd);

out = struct('yff', zeros(N.t,1), 'yft', zeros(M,1),...
             'ytf', zeros(M,1), 'ytt', zeros(N.t,1),...
             'ysh', zeros(N.t,1), 'S', zeros(N.t,1), 'pqopt', pqopt);
if vtau
    out.lntau = log(tau);
    out.taupq = exp(u0(bidx.pq)).*cnx.F.pq'*(yff.*tau.*log(tau));
    out.taupv = exp(u0(bidx.pv)).*cnx.F.pv'*(yff.*tau.*log(tau));
    out.delta = tshift;
%     out.delta = zeros(size(tshift));
else
    out.lntau = zeros(M,1);
    out.taupq = zeros(N.pq,1);
    out.taupv = zeros(N.pv,1);
    out.delta = zeros(size(tshift));
end
%% PV buses
switch pqopt
    case 1
        %% PQ bues (option 1)
        out.yff(bidx.pv)    = exp(2*u0(bidx.pv)).*cnx.F.pv'*(yff.*tau);
        out.ytt(bidx.pv)    = exp(2*u0(bidx.pv)).*cnx.T.pv'*ytt;

        ytmp                = yft.*exp((F+T)*u0).*exp(-1i*tshift);
        out.yft(lidx.F.pv)  = ytmp(lidx.F.pv);

        ytmp                = ytf.*exp((T+F)*u0).*tau.*exp(1i*tshift);
        out.ytf(lidx.T.pv)  = ytmp(lidx.T.pv);

    case 2
        %% PQ buses (option 2)
        out.yff(bidx.pv)    = exp(2*u0(bidx.pv)).*cnx.F.pv'*yff;
        out.ytt(bidx.pv)    = exp(2*u0(bidx.pv)).*cnx.T.pv'*ytt;

        ytmp                = exp((F+T)*u0).*(yft./tau).*exp(-1i*tshift);
        out.yft(lidx.F.pv)  = ytmp(lidx.F.pv);

        ytmp                = exp((F+T)*u0).*ytf.*tau.*exp(1i*tshift);
        out.ytf(lidx.T.pv)  = ytmp(lidx.T.pv);
        
    otherwise
        error('reweightY: pqopt must be either 1 or 2')
end
% out.yff(bidx.pv)    = exp(2*u0(bidx.pv)).*cnx.F.pv'*yff;
% out.ytt(bidx.pv)    = exp(2*u0(bidx.pv)).*cnx.T.pv'*ytt;
% 
% ytmp                = yft.*exp((F+T)*u0).*exp(-1i*tshift)./tau;
% % ytmp(lidx.T.pq)     = ytmp(lidx.T.pq)./tau(lidx.T.pq); % vf/tau option
% out.yft(lidx.F.pv)  = ytmp(lidx.F.pv);
% 
% ytmp                = ytf.*exp((F+T)*u0).*exp(1i*tshift).*tau;
% % ytmp(lidx.F.pq)     = ytmp(lidx.F.pq).*tau(lidx.F.pq); % vf/tau option
% out.ytf(lidx.T.pv)  = ytmp(lidx.T.pv);
% 
out.ysh(bidx.pv)    = exp(2*u0(bidx.pv)).*ysh(bidx.pv);
out.S(bidx.pv)      = Sp.Pg(bidx.pv) - Sp.Pd(bidx.pv) - 1i*Sp.Qd(bidx.pv);
%% PQ buses
switch pqopt
    case 1
        %% PQ bues (option 1)
        out.yff(bidx.pq)    = exp(u0(bidx.pq)).*cnx.F.pq'*(yff.*tau);
        out.ytt(bidx.pq)    = exp(u0(bidx.pq)).*cnx.T.pq'*ytt;

        ytmp                = yft.*exp(T*u0).*exp(-1i*tshift);
        out.yft(lidx.F.pq)  = ytmp(lidx.F.pq);

        ytmp                = ytf.*exp(F*u0).*tau.*exp(1i*tshift);
        out.ytf(lidx.T.pq)  = ytmp(lidx.T.pq);

        out.ysh(bidx.pq)    = exp(u0(bidx.pq)).*ysh(bidx.pq);
        out.S(bidx.pq)      = exp(-u0(bidx.pq)).*S(bidx.pq);
    case 2
        %% PQ buses (option 2)
        out.yff(bidx.pq)    = cnx.F.pq'*yff;
        out.ytt(bidx.pq)    = cnx.T.pq'*ytt;

        ytmp                = exp(-(F-T)*u0).*(yft./tau).*exp(-1i*tshift);
        out.yft(lidx.F.pq)  = ytmp(lidx.F.pq);

        ytmp                = exp((F-T)*u0).*ytf.*tau.*exp(1i*tshift);
        out.ytf(lidx.T.pq)  = ytmp(lidx.T.pq);

        out.ysh(bidx.pq)    = ysh(bidx.pq);
        out.S(bidx.pq)      = exp(-2*u0(bidx.pq)).*S(bidx.pq);
    otherwise
        error('reweightY: pqopt must be either 1 or 2')
end