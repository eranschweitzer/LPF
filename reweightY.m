function out = reweightY(opmats,Sb,Sp,u0,varargin)

pqopt = varargin_parse(varargin,'pqopt',1);
v2struct(opmats); %move fields to variables
%% define base values
ys = Sb.g + 1i*Sb.b;
ytt = ys + 0.5*1i*Sb.bc;
yff = (ytt)./Sb.tau.^2;
yft = -ys.*(exp(+1i*Sb.tshift)./Sb.tau);
ytf = -ys.*(exp(-1i*Sb.tshift)./Sb.tau);
ysh = Sb.gsh + 1i*Sb.bsh;
S   = (Sp.Pg - Sp.Pd) + 1i*(Sp.Qg - Sp.Qd);

out = struct('yff', zeros(N.t,1), 'yft', zeros(M,1),...
             'ytf', zeros(M,1), 'ytt', zeros(N.t,1),...
             'ysh', zeros(N.t,1), 'S', zeros(N.t,1), 'pqopt', pqopt);
%% PV buses
out.yff(bidx.pv)    = exp(2*u0(bidx.pv)).*cnx.F.pv'*yff;
out.ytt(bidx.pv)    = exp(2*u0(bidx.pv)).*cnx.T.pv'*ytt;

ytmp                = yft.*exp((F+T)*u0);
out.yft(lidx.F.pv)  = ytmp(lidx.F.pv);

ytmp                = ytf.*exp((F+T)*u0);
out.ytf(lidx.T.pv)  = ytmp(lidx.T.pv);

out.ysh(bidx.pv)    = exp(2*u0(bidx.pv)).*ysh(bidx.pv);
out.S(bidx.pv)      = Sp.Pg(bidx.pv) - Sp.Pd(bidx.pv) - 1i*Sp.Qd(bidx.pv);

switch pqopt
    case 1
        %% PQ bues (option 1)
        out.yff(bidx.pq)    = exp(u0(bidx.pq)).*cnx.F.pq'*yff;
        out.ytt(bidx.pq)    = exp(u0(bidx.pq)).*cnx.T.pq'*ytt;

        ytmp                = yft.*exp(T*u0);
        out.yft(lidx.F.pq)  = ytmp(lidx.F.pq);

        ytmp                = ytf.*exp(F*u0);
        out.ytf(lidx.T.pq)  = ytmp(lidx.T.pq);

        out.ysh(bidx.pq)    = exp(u0(bidx.pq)).*ysh(bidx.pq);
        out.S(bidx.pq)      = exp(-u0(bidx.pq)).*S(bidx.pq);
    case 2
        %% PQ buses (option 2)
        out.yff(bidx.pq)    = cnx.F.pq'*yff;
        out.ytt(bidx.pq)    = cnx.T.pq'*ytt;

        ytmp                = yft.*exp(-(F-T)*u0);
        out.yft(lidx.F.pq)  = ytmp(lidx.F.pq);

        ytmp                = ytf.*exp((F-T)*u0);
        out.ytf(lidx.T.pq)  = ytmp(lidx.T.pq);

        out.ysh(bidx.pq)    = ysh(bidx.pq);
        out.S(bidx.pq)      = exp(-2*u0(bidx.pq)).*S(bidx.pq);
    otherwise
        error('reweightY: pqopt must be either 1 or 2')
end