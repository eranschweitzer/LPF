function vars = solve_iterphi(ids,F,T,E,Sb,Sp,u0,bidx,theta_ref,N,varargin)

epsilon = varargin_parse(varargin,'epsilon',1e-8);
itermax = varargin_parse(varargin,'itermax',25);
% GwAru = varargin_parse(varargin,'GwAru','default');
% GwAmu = varargin_parse(varargin,'GwAmu','default');
Gw = varargin_parse(varargin,'Gw',branchweights(Sb));
phi = varargin_parse(varargin,'phi',0);

%% matrices fixing PV and Ref quantities
[I,U] = pv_ref_mats(bidx,N,u0);

%% INITIALIZATION
% phi = zeros(size(E,1),1);
[Art,Aru,Amt,Amu,br,bm] = matrix_parts_init(ids,F,T,E,u0,Sb,Sp,'Gw',Gw,'phi',phi);

%% solution
maxerr = struct('phi',0);
convg = 0;
for k = 0:itermax    
    x      = single_solve(Art,Aru,Amt,Amu,br,bm,I,U,theta_ref,N);
    vars   = result_parse(x,u0,N);
    phinew = 0.5*(E*vars.theta).^2;
    maxerr.phi = max(abs(phinew - phi));
    if any(isnan(phinew)) || any(isinf(phinew))
        % solution diverged
        break
    elseif maxerr.phi <= epsilon
        convg = 1;
        break
    else
        phi = phinew;
        % update br and bm
        br  = genbr(ids.br,ids.brb,F,T,Sb,Sp.Pg,Sp.Pd,u0,phi,'Gw',Gw);
        bm  = genbm(ids.bm,ids.bmb,F,T,Sb,Sp.Qg,Sp.Qd,u0,phi,'Gw',Gw); 
    end
end

vars.phi = phinew;
vars.v  = exp(vars.u0 + vars.uhat);
vars.convg = convg;
vars.iter = k;
vars.maxerr = maxerr;
