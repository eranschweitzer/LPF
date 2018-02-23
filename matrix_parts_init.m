function [Art,Aru,Amt,Amu,br,bm] = matrix_parts_init(ids,F,T,E,u0,Sb,Sp,varargin)

% GwAru = varargin_parse(varargin,'GwAru','default');
% GwAmu = varargin_parse(varargin,'GwAmu','default');
Gw  = varargin_parse(varargin,'Gw',branchweights(Sb));
phi = varargin_parse(varargin,'phi',0);

if length(fieldnames(ids)) == 13
    Art = genArtheta(ids.Art,ids.Artb,F,T,E,Sb,u0);
    Aru = genAru(ids.Aru,ids.Arub,F,T,Sb,Sp.Pg,Sp.Pd,u0, 'Gw',Gw);
    Amt = genAmtheta(ids.Amt,ids.Amtb,F,T,E,Sb,u0);
    Amu = genAmu(ids.Amu,ids.Amub,F,T,Sb,Sp.Qg,Sp.Qd,u0, 'Gw',Gw);
    br  = genbr(ids.br,ids.brb,F,T,Sb,Sp.Pg,Sp.Pd,u0,phi,'Gw',Gw);
    bm  = genbm(ids.bm,ids.bmb,F,T,Sb,Sp.Qg,Sp.Qd,u0,phi,'Gw',Gw);
elseif length(fieldnames(ids)) == 12
    betas = genBetas(ids,F,T,E,u0,Sb,Sp);
    [lambda, psi, nu] = genUcoeffs(ids, F, T, E, u0, Sb, Sp);
    gammas = genGammas(ids, F, T, E, u0, Sb, Sp);
    
    br  = betas.sh.r - F'*betas.f.r - T'*betas.t.r;
    bm  = betas.sh.m - F'*betas.f.m - T'*betas.t.m;
    Art = -(F'*diags(gammas.f.m)*E + T'*diags(gammas.t.m)*E);
    Aru = diags(lambda.r) + F'*(diags(psi.f.r)*F - diags(psi.t.r)*T) + T'*(diags(nu.f.r)*F - diags(nu.t.r)*T);
    Amt = F'*diags(gammas.f.r)*E + T'*diags(gammas.t.r)*E;
    Amu = diags(lambda.m) + F'*(diags(psi.f.m)*F - diags(psi.t.m)*T) + T'*(diags(nu.f.m)*F - diags(nu.t.m)*T);
    
    if ids.c < 2
        br  = br + Aru*u0;
        bm  = bm + Amu*u0;
    end
else
    error('Incorrect Ids structures. Number of fields is %d',  length(fieldnames(ids)))
end
