function vars = solve_iterphiu0(ids,F,T,E,Sb,Sp,u0,bidx,theta_ref,N,varargin)

epsilon_phi = varargin_parse(varargin,'epsilon_phi',1e-8);
epsilon_u0  = varargin_parse(varargin,'epsilon_u0',1e-8);
itermax = varargin_parse(varargin,'itermax',25);
% GwAru = varargin_parse(varargin,'GwAru','default');
% GwAmu = varargin_parse(varargin,'GwAmu','default');
Gw = varargin_parse(varargin,'Gw',branchweights(Sb));
phi = varargin_parse(varargin,'phi',0);
uopt = varargin_parse(varargin,'uopt','flat');
%% matrices fixing PV and Ref quantities
[I,U] = pv_ref_mats(bidx,N,u0);
%% INITIALIZATION
nv = [2,4,5,7]; %cases (IDs) where u0 is not considered
% phi = zeros(size(E,1),1);
[Art,Aru,Amt,Amu,br,bm] = matrix_parts_init(ids,F,T,E,u0,Sb,Sp,'Gw',Gw,'phi',phi);

%% solution
maxerr = struct('phi',0,'u0',0);
convg = 0;
for k = 0:itermax    
    x      = single_solve(Art,Aru,Amt,Amu,br,bm,I,U,theta_ref,N);
    vars   = result_parse(x,u0,N);
    phinew = 0.5*(E*vars.theta).^2;
    u0new  = vars.uhat + u0;
    u0new  = vlimit(u0new,E);
    if k < 4
        u0new  = vlimit(u0new,E,'uopt',uopt);
    else
        % fix nodes to pv buses
        [u0new,busmask]  = vlimit(u0new,E,'uopt',uopt);
        if any(busmask)
            [bidx,N,I,U] = fixbus(busmask,u0new,bidx,N);
        end
    end
    maxerr.phi = max(abs(phinew - phi));
    maxerr.u0  = max(abs(u0new - u0));
    if any(isnan(u0new)) || any(isinf(u0new)) || any(isnan(phinew)) || any(isinf(phinew))
        % solution diverged
        break
    elseif (maxerr.phi <= epsilon_phi) && (maxerr.u0 <= epsilon_u0)
        convg = 1;
        break
    else
        phi = phinew;
        u0 = u0new;
        % update br and bm
        br  = genbr(ids.br,ids.brb,F,T,Sb,Sp.Pg,Sp.Pd,u0,phi,'Gw',Gw);
        bm  = genbm(ids.bm,ids.bmb,F,T,Sb,Sp.Qg,Sp.Qd,u0,phi,'Gw',Gw); 
        % update parts involving u0
        if ~ismember(ids.Art,nv)
            Art = genArtheta(ids.Art,ids.Artb,F,T,E,Sb,u0);
        end
        if ~ismember(ids.Aru,nv)
            Aru = genAru(ids.Aru,ids.Arub,F,T,Sb,Sp.Pg,Sp.Pd,u0,'Gw',Gw);
        end
        if ~ismember(ids.Amt,nv)
            Amt = genAmtheta(ids.Amt,ids.Amtb,F,T,E,Sb,u0);
        end
        if ~ismember(ids.Amu,nv)
            Amu = genAmu(ids.Amu,ids.Amub,F,T,Sb,Sp.Qg,Sp.Qd,u0,'Gw',Gw);
        end
    end
end

vars.phi  = phinew;
vars.u0   = u0new;
vars.uhat = u0new - u0;
vars.v  = exp(vars.u0 + vars.uhat);
vars.convg = convg;
vars.iter = k;
vars.maxerr = maxerr;
 