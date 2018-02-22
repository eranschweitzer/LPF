function [lambda, psi, nu] = genUcoeffs(ids, F, T, E, u0, Sb, Sp, varargin)

branch_flag = varargin_parse(varargin, 'branch', false);
only   = varargin_parse(varargin, 'only', '');

if ~branch_flag
	Pg = Sp.Pg; Pd = Sp.Pd; 
	Qg = Sp.Qg; Qd = Sp.Qd;
	if ids.c == 2
		ctmp0 = ones(size(u0));
		ctmp = ones(size(u0));
	elseif ids.c == 1
		ctmp0 = exp(u0);
		ctmp = ones(size(u0));
	else
		ctmp0 = exp(u0);
		ctmp = exp(u0);
	end
else
	Pg = 0; Pd = 0; 
	Qg = 0; Qd = 0;
	if ids.c > 0
		ctmp0 = ones(size(u0));
		ctmp  = ones(size(u0));
	else
		ctmp0 = exp(u0);
		ctmp  = exp(u0);
	end
end

%% real
if strcmp(only,'') || strcmp(only,'real')
	if ~branch_flag
		[g,b,~,tau,tshift,gsh,~] = branchparts_unpack(Sb,ids.gr>0);
		if ids.Omr == 1
			% simplification means setting tau = 1
			tau = ones(size(tau));
		end
		
		if ids.gr > 0
			g = zeros(size(g));
		end
		% multiplying by decoup will zero out coefficients if ids.gr = 2
		% thus decoupling real(K) from u
		if ids.gr == 2
			decoup = 0;
		else
			decoup = 1;
		end
	else
		[g,b,~,tau,tshift,gsh,~] = branchparts_unpack(Sb,ids.g);
		gsh = 0;
		if ids.g == 1
			g = zeros(size(g));
			decoup = 0;
		else
			decoup = 1;
		end
		if ids.simp == 1
			% simplification means setting tau = 1
			tau = ones(size(tau));
		end
	end
	
	%%%%%%% lambda %%%%%%%%%%%
	if ~branch_flag
		lambda.r = (Pg - Pd).*ctmp0.^(-1) + gsh.*ctmp0;
		lambda.r = lambda.r*decoup;
	else
		lambda.r = 0;
	end
	
	%%%%%%% psi and nu  %%%%%%%%%%%%%
	
	if ids.AB == 1
		psi.f.r = g./tau.*(F*ctmp)./tau;
		psi.t.r = g./tau.*(T*ctmp).*cos(tshift) - b./tau.*(T*ctmp).*sin(tshift);
		psi.t.r = psi.t.r*decoup;
		
		nu.f.r  = -g./tau.*(F*ctmp).*cos(tshift) - b./tau.*(F*ctmp).*sin(tshift);
		nu.f.r  = nu.f.r*decoup;
		nu.t.r  = -g.*(T*ctmp);
	else
		psi.f.r = g./tau.*(F*ctmp);
		psi.t.r = g./tau.*(T*ctmp);
		nu.f.r = -g.*(F*ctmp);
		nu.t.r = -g.*(T*ctmp);
	end
end

%% imaginary
if strcmp(only,'') || strcmp(only,'imag')
	if ~branch_flag
		[g,b,bc,tau,tshift,~,bsh] = branchparts_unpack(Sb,ids.gi>0);
		if ids.Omi == 1
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
		%%%%%% labmda %%%%%%%%%%%
		if ~branch_flag
			lambda.m = (Qd - Qg).*ctmp0.^(-1) + ctmp0.*(bsh + F'*(b0./tau.^2) + T'*b0);
		else
			lambda.m.f = (F*ctmp0).*b0./tau.^2;
			lambda.m.t = b0;
		end
		%%%%%%%%%%% psi and nu %%%%%%%%%%%
		if ids.AB == 1
			psi.f.m = b./tau.*(F*ctmp)./tau;
			psi.t.m = b./tau.*(T*ctmp).*cos(tshift) + g./tau.*(T*ctmp).*sin(tshift);
	
			nu.f.m = -b./tau.*(F*ctmp).*cos(tshift) + g./tau.*(F*ctmp).*sin(tshift);
			nu.t.m = -b.*(T*ctmp);		
		else
			psi.f.m = b./tau.*(F*ctmp);
			psi.t.m = b./tau.*(T*ctmp);
	
			nu.f.m = -b.*(F*ctmp);
			nu.t.m = -b.*(T*ctmp);
		end
	else
		%%%%%% labmda %%%%%%%%%%%
		if ~branch_flag
			lambda.m = 2*(Qd - Qg).*ctmp0.^(-2); 
		else
			lambda.m.f = 0; lambda.m.t = 0;
		end
		%%%%%%%%%%% psi and nu %%%%%%%%%%%
		if ids.AB == 1
			psi.f.m = b./tau.*cos(tshift) + g./tau.*sin(tshift);
			psi.t.m = psi.f.m;
	
			nu.f.m = -b./tau.*cos(tshift) + g./tau.*sin(tshift);
			nu.t.m = nu.f.m;
		else
			psi.f.m = b./tau.^(2);
			psi.t.m = psi.f.m;
	
			nu.f.m  = -b;
			nu.t.m  = nu.f.m;
		end
	end
end
