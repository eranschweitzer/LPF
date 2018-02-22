function mat = genAmtheta(ID,btype,F,T,E,Sb,u0)

[g,b,~,tau,tshift,~,~] = branchparts_unpack(Sb,btype);
N = length(u0);
switch ID
    case 0
        mat = ...
        (F'*diags(tau.^(-1))*diags(exp(T*u0).*(g.*cos(tshift) - b.*sin(tshift)))...
        -T'*diags(tau.^(-1))*diags(exp(F*u0).*(g.*cos(tshift) + b.*sin(tshift))))*E;
    case 1
        mat = ...
        (F'*diags(tau.^(-1))*diags(exp(T*u0).*g)...
        -T'*diags(tau.^(-1))*diags(exp(F*u0).*g))*E;
    case 2
        mat = ...
        (F'*diags(tau.^(-1))*diags(g)...
        -T'*diags(tau.^(-1))*diags(g))*E;
    case 5
        mat = sparse(N,N);
    otherwise
        error('Incorrect ID')
end