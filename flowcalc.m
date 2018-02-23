function S = flowcalc(id, pq, vars, F, T, E, Sb)

if id.Q < 3
	% if id.lin = 1 (1-id.lin)=0 and both af and at will be simply a vector of 1s linearizing the problem
	% the factor id.Q corresponds exactly to wether we have exp(1u) or exp(2u)
  % note that u0 + uhat = u	
	af = exp((1-id.lin)*(id.Q)*F*vars.u);
	at = exp((1-id.lin)*(id.Q)*T*vars.u);

	betas = genBetas(id,F,T,E,vars.u0,Sb,0, 'branch', true, 'only', pq);
	[lambda, psi, nu] = genUcoeffs(id, F, T, E, vars.u0, Sb, 0, 'branch', true, 'only', pq);
	gammas = genGammas(id, F, T, E, vars.u0, Sb, 0, 'branch', true, 'only', pq);
	
    if id.c == 2
        u0tmp.self = zeros(size(vars.u0));
        u0tmp.cross= zeros(size(vars.u0));
    elseif id.c == 1
        u0tmp.self = zeros(size(vars.u0));
        u0tmp.cross= vars.u0;
    elseif id.c == 0
        u0tmp.self = vars.u0;
        u0tmp.cross= vars.u0;
    end
    
	S = struct();
	if strcmp(pq,'real')
        S.f = af.*(betas.f.r + (diags(psi.f.r)*F - diags(psi.t.r)*T)*(vars.u - u0tmp.cross) - gammas.f.m.*(E*vars.theta) );
		S.t = at.*(betas.t.r + (diags(nu.f.r)*F  - diags(nu.t.r)*T) *(vars.u - u0tmp.cross) - gammas.t.m.*(E*vars.theta) );
% 		S.f = af.*(betas.f.r + psi.f.r.*(F*vars.u) - psi.t.r.*(T*vars.u) - gammas.f.m.*(E*vars.theta) );
% 		S.t = at.*(betas.t.r - nu.t.r.*(T*vars.u)  + nu.f.r.*(F*vars.u)  - gammas.t.m.*(E*vars.theta) );
	elseif strcmp(pq, 'imag')
        S.f = -af.*( betas.f.m + diags(lambda.m.f)*F*(vars.u - u0tmp.self) + ...
            (diags(psi.f.m)*F - diags(psi.t.m)*T)*(vars.u - u0tmp.cross) + gammas.f.r.*(E*vars.theta) );
		S.t = -at.*( betas.t.m + diags(lambda.m.t)*T*(vars.u - u0tmp.self) + ...
            (diags(nu.f.m)*F  - diags(nu.t.m)*T)*(vars.u - u0tmp.cross)  + gammas.t.r.*(E*vars.theta) );
% 		S.f = -af.*( betas.f.m + (lambda.m.f + psi.f.m).*(F*vars.u) - psi.t.m.*(T*vars.u) + gammas.f.r.*(E*vars.theta) );
% 		S.t = -at.*( betas.t.m + (lambda.m.t - nu.t.m).*(T*vars.u)  + nu.f.m.*(F*vars.u)  + gammas.t.r.*(E*vars.theta) );
	else
		error('pq must be either real or imag %s given', pq)
	end
elseif id.Q == 3
	% Zhigang's calculation
	if strcmp(pq, 'imag') || strcmp(pq, 'real')
		S = zhigangflow(vars, F, T, E, Sb, 'pq', pq);
	else
		error('pq must be either real or imag %s given', pq)
	end
elseif id.Q == 4
	% full nonlinear AC
	if strcmp(pq, 'real')
		S = calcPflow(5, 0, vars, F, T, E, Sb);
	elseif strcmp(pq, 'imag')
		S = calcQflow(5, 0, vars, F, T, E, Sb);
	else
		error('pq must be either real or imag %s given', pq)
	end
end
