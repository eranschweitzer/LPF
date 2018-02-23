function betas = genBetas(ids, F, T, E, u0, Sb, Sp, varargin)

branch_flag = varargin_parse(varargin, 'branch', false);
only   = varargin_parse(varargin, 'only', '');

if ~branch_flag
	Pg = Sp.Pg; Pd = Sp.Pd; 
	Qg = Sp.Qg; Qd = Sp.Qd;
	if ids.c == 2
		ctmp0 = ones(size(u0));
		ctmp  = ones(size(u0));
		%u0tmp = zeros(size(u0));
	elseif ids.c == 1
		ctmp0 = exp(u0);
		ctmp = ones(size(u0));
		%u0tmp = zeros(size(u0));
	else
		ctmp0 = exp(u0);
		ctmp  = exp(u0);
		%u0tmp = u0;
	end
else
	Pg = 0; Pd = 0; 
	Qg = 0; Qd = 0;
	if ids.c > 0
		ctmp0 = ones(size(u0));
		ctmp  = ones(size(u0));
		%u0tmp = zeros(size(u0));
	else
		ctmp0 = exp(u0);
		ctmp  = exp(u0);
		%u0tmp = u0;
	end
end

%% real part
if strcmp(only,'') || strcmp(only,'real')
	if ~branch_flag
		[g,b,~,tau,tshift,gsh,~] = branchparts_unpack(Sb,ids.gr>0);
		if ids.Kr == 1
			% simplification means setting tau = 1
			tau = ones(size(tau));
		end
		
		if ids.gr > 0
			g = zeros(size(g));
		end
	else
		[g,b,~,tau,tshift,gsh,~] = branchparts_unpack(Sb,ids.g);
		gsh = 0;
		if ids.g == 1
			g = zeros(size(g));
		end
		if ids.simp == 1
			% simplification means setting tau = 1
			tau = ones(size(tau));
		end
	end
	%%%% self %%%%%%
	if ~branch_flag
		betas.sh.r = (Pg - Pd).*ctmp0.^(-1) - gsh.*ctmp;
	end
	%%%% cross %%%%%%
	if ids.AB == 1
		betas.f.r = -g./tau.*(-(F*ctmp)./tau + (T*ctmp).*cos(tshift)) + b./tau.*(T*ctmp).*sin(tshift);
		betas.t.r = g.*( (T*ctmp) - (F*ctmp).*cos(tshift)./tau) - b./tau.*(F*ctmp).*sin(tshift);
	else
		betas.f.r = -g./tau.*( (F*ctmp).*log(tau) - (E*ctmp)) + b./tau.*(T*ctmp).*tshift;
		betas.t.r = g.*( (F*ctmp).*log(tau) - (E*ctmp) ) - b.*(F*ctmp).*tshift;
	end
end

%% imaginary part
if strcmp(only,'') || strcmp(only,'imag')
	if ~branch_flag
		[g,b,bc,tau,tshift,~,bsh] = branchparts_unpack(Sb,ids.gi>0);
		if ids.Ki == 1
			% simplification means setting tshift = 0
			tshift = zeros(size(tshift));
		end
		if ids.gi == 1
			g = zeros(size(g));
		end
	else
		[g,b,bc,tau,tshift,~,bsh] = branchparts_unpack(Sb,ids.g);
		bsh = 0;
		if ids.g == 0
			g = zeros(size(g));
		end
		if ids.simp == 1
			% simplification means setting tshift = 0
			tshift = zeros(size(tshift));
		end
	end
	b0 = bc/2;
	
	if ids.Q == 1
		%%%%%%%% self %%%%%%%%
		if ~branch_flag
			betas.sh.m = (Qd - Qg).*ctmp0.^(-1) - bsh.*ctmp0;
		end
		%%%%%%% cross %%%%%%%%
		if ids.AB == 1
				betas.f.m = (b0./(tau.^2)).*(F*ctmp) - b./tau.*( -(F*ctmp)./tau + (T*ctmp).*cos(tshift) ) - g./tau.*(T*ctmp).*sin(tshift);
				betas.t.m = b0.*(T*ctmp) + b.*( (T*ctmp) - (F*ctmp).*cos(tshift)./tau) + g./tau.*(F*ctmp).*sin(tshift);
		else
				betas.f.m = b0./tau.*(F*ctmp).*(1 - log(tau)) - b./tau.*( (F*ctmp).*log(tau) - (E*ctmp) ) - g./tau.*(T*ctmp).*tshift;
				betas.t.m = b0.*(T*ctmp) + b.*( (F*ctmp).*log(tau) - (E*ctmp) ) + g.*(F*ctmp).*tshift;
		end
	else
		%%%%%%% self %%%%%%%%%
		betas.sh.m = (Qd - Qg).*ctmp0.^(-2) - bsh;
		%%%%%% cross %%%%%%%
		ctf = (T*ctmp)./(F*ctmp);
		cft = (F*ctmp)./(T*ctmp);
		if ids.AB == 1
				betas.f.m = (b0./(tau.^2)) - b./(tau.^2).*( ctf.*tau.*cos(tshift) - 1) - g./tau.*sin(tshift).*ctf;
				betas.t.m = b0 + b.*(1 - cft.*cos(tshift)./tau ) + g./tau.*sin(tshift).*cft;
				%betas.f.m = (b0./(tau.^2)) - b./(tau.^2).*(tau.*cos(tshift).*(1 - (E*u0tmp)) - 1) - g./tau.*sin(tshift).*(1 - (E*u0tmp));
				%betas.t.m = b0 + b.*(1 - cos(tshift)./tau.*(1 + (E*u0tmp) ) ) + g./tau.*sin(tshift).*(1 + (E*u0tmp) );
		else
				betas.f.m = b0./(tau.^2) - b./(tau.^2).*(ctf.*(1 + log(tau)) - 1 ) - g./(tau.^2).*tshift.*ctf;
				betas.t.m = b0 + b.*( 1 - cft.*( 1 - log(tau)) ) + g.*tshift.*cft;
				%betas.f.m = b0./(tau.^2) - b./(tau.^2).*(log(tau) - (E*u0tmp) ) - g./(tau.^2).*tshift;
				%betas.t.m = b0 + b.*( log(tau) - (E*u0tmp) ) + g.*tshift;
		end
	end
end
