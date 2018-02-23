function vars = pfsolve(ids,F,T,E,Sb,Sp,bidx,theta_ref,N,u0,varargin)

%% optional inputs
epsilon = varargin_parse(varargin,'epsilon',1e-8);
itermax = varargin_parse(varargin,'itermax',25);
% GwAru = varargin_parse(varargin,'GwAru','default');
% GwAmu = varargin_parse(varargin,'GwAmu','default');
Gw  = varargin_parse(varargin,'Gw',branchweights(Sb));
phi = varargin_parse(varargin,'phi',0);
if ids.iter == 0
    uopt = varargin_parse(varargin,'uopt','none');
else
    uopt = varargin_parse(varargin,'uopt','flat');
end
%%
switch ids.iter
    case 0
        % No iteration
        [Art,Aru,Amt,Amu,br,bm] = matrix_parts_init(ids,F,T,E,u0,Sb,Sp,'Gw',Gw,'phi',phi);
        [matI,matU] = pv_ref_mats(bidx,N,u0,ids);
        x    = single_solve(Art,Aru,Amt,Amu,br,bm,matI,matU,theta_ref,N);
        if length(fieldnames(ids)) == 12 
            vars = result_parse(x,u0,N,'uvar', 'u');
            vars.v     = exp(vars.u);
        else
            vars = result_parse(x,u0,N);
            u    = vlimit(vars.u0+vars.uhat,E,'uopt',uopt);
            vars.uhat  = u - vars.u0;
            vars.v     = exp(u);
        end
        vars.convg = 1;      
        vars.phi   = phi;
    case 1
        % Only phi iteration
        vars = solve_iterphi(ids,F,T,E,Sb,Sp,u0,bidx,theta_ref,N,...
            'epsilon',epsilon,'itermax',itermax,'Gw',Gw,'phi',phi,'uopt',uopt);
    case 2
        % Only u0 iteration
        vars = solve_iteru0(ids,F,T,E,Sb,Sp,u0,bidx,theta_ref,N,...
            'epsilon',epsilon,'itermax',itermax,'Gw',Gw,'phi',phi,'uopt',uopt);
    case 3
        % Phi and u0 iteration
        vars = solve_iterphiu0(ids,F,T,E,Sb,Sp,u0,bidx,theta_ref,N,...
            'epsilon',epsilon,'itermax',itermax,'Gw',Gw,'phi',phi,'uopt',uopt);
    otherwise
        error('Incorrect Iterration ID')
end
