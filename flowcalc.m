function S = flowcalc(id, pq, vars, F, T, E, Sb)

if id.Q < 3
	% if id.lin = 1 (1-id.lin)=0 and both af and at will be simply a vector of 1s linearizing the problem
	% the factor id.Q corresponds exactly to wether we have exp(1u) or exp(2u)
  % note that u0 + uhat = u	
	af = exp((1-id.lin)*(id.Q)*F*(vars.u0 + vars.uhat));
	at = exp((1-id.lin)*(id.Q)*T*(vars.u0 + vars.uhat));

	betas = genBetas(id,F,T,E,vars.u0,Sb,0, 'branch', true, 'only', pq);
	[lambda, psi, nu] = genUcoeffs(id, F, T, E, vars.u0, Sb, 0, 'branch', true, 'only', pq);
	gammas = genGammas(id, F, T, E, vars.u0, Sb, 0, 'branch', true, 'only', pq);
	
	S = struct();
	if strcmp(pq,'real')
		S.f = af.*(betas.f.r + psi.f.r.*(F*vars.uhat) - psi.t.r.*(T*vars.uhat) - gammas.f.m.*(E*vars.theta) );
		S.t = at.*(betas.t.r - nu.t.r.*(T*vars.uhat)  + nu.f.r.*(F*vars.uhat)  - gammas.t.m.*(E*vars.theta) );
	elseif strcmp(pq, 'imag')
		S.f = -af.*( betas.f.m + (lambda.m.f + psi.f.m).*(F*vars.uhat) - psi.t.m.*(T*vars.uhat) + gammas.f.r.*(E*vars.theta) );
		S.t = -at.*( betas.t.m + (lambda.m.t - nu.t.m).*(T*vars.uhat)  + nu.f.m.*(F*vars.uhat)  + gammas.t.r.*(E*vars.theta) );
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
