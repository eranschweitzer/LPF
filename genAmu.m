function mat = genAmu(ID,btype,F,T,Sb,Qg,Qd,u0,varargin)

[g,b,bc,tau,tshift,~,bsh] = branchparts_unpack(Sb,btype);
% Gwtype = varargin_parse(varargin,'Gw','default');
% Gw = genweights(b,Gwtype);
Gw = varargin_parse(varargin,'Gw',branchweights(Sb));
switch ID
    case 0
        mat = ...
        diags(exp(-u0).*(Qd-Qg) + exp(u0).*bsh)...
        +F'*diags(tau.^(-1))*( diags(tau.^(-1))*diags(exp(Gw.b*F*u0).*(b+bc/2))*F...
           -diags(exp(Gw.b*T*u0).*(g.*sin(tshift) + b.*cos(tshift)))*T)...
        +T'*(diags(exp(Gw.b*T*u0).*(b+bc/2))*T ...
           -diags(tau.^(-1))*diags(exp(Gw.b*F*u0).*(-g.*sin(tshift) + b.*cos(tshift)))*F);
    case 1
        mat = ...
        diags(exp(-u0).*(Qd-Qg) + exp(u0).*bsh)...
        +F'*diags(tau.^(-1))*( diags(tau.^(-1))*diags(exp(Gw.b*F*u0).*(b+bc/2))*F...
           -diags(exp(Gw.b*T*u0).*b)*T)...
        +T'*(diags(exp(Gw.b*T*u0).*(b+bc/2))*T ...
           -diags(tau.^(-1))*diags(exp(Gw.b*F*u0).*b)*F);
    case 2
        mat = ...
        diags(Qd - Qg + bsh)...
        +F'*diags(tau.^(-1))*( diags(tau.^(-1))*diags(b+bc/2)*F - diags(b)*T)...
        +T'*(diags(b+bc/2)*T -diags(tau.^(-1))*diags(b)*F);
    otherwise
        error('Incorrect ID')
end