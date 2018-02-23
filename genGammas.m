function gammas = genGammas(ids, F, T, E, u0, Sb, Sp, varargin)

branch_flag = varargin_parse(varargin, 'branch', false);
only   = varargin_parse(varargin, 'only', '');

if ids.c > 0
	ctmp = ones(size(u0));
else
	ctmp = exp(u0);
end
%% real
if strcmp(only,'') || strcmp(only,'imag')
	% Recall that the real part of Gamma affects the imaginary parts of K
	
	if ~branch_flag
		[g,b,~,tau,tshift,gsh,~] = branchparts_unpack(Sb,ids.gi > 0);
		if ids.Gamr == 1
			% simplification means setting tshift = 0
			tshift = zeros(size(tshift));
		end
		if ids.gi > 0
			g = zeros(size(g));
		end
	else
		[g,b,~,tau,tshift,gsh,~] = branchparts_unpack(Sb,ids.g);
		if ids.g > 0
			g = zeros(size(g));
		end
		if ids.simp == 1
			% simplification means setting tshift = 0
			tshift = zeros(size(tshift));
		end
	end
	
	if ids.Q == 1
		if ids.AB == 1
			gammas.f.r =  g./tau.*(T*ctmp).*cos(tshift) - b./tau.*(T*ctmp).*sin(tshift);
			gammas.t.r = -g./tau.*(F*ctmp).*cos(tshift) - b./tau.*(F*ctmp).*sin(tshift);
		else
			gammas.f.r =  g./tau.*(T*ctmp);
			gammas.t.r = -g.*(F*ctmp);
		end
	else
		ctf = (T*ctmp)./(F*ctmp);
		cft = (F*ctmp)./(T*ctmp);
		if ids.AB == 1
			gammas.f.r =  ctf.*( g./tau.*cos(tshift) - b./tau.*sin(tshift) );
			gammas.t.r = -cft.*( g./tau.*cos(tshift) + b./tau.*sin(tshift) );
			%gammas.f.r =  g./tau.*cos(tshift) - b./tau.*sin(tshift);
			%gammas.t.r = -g./tau.*cos(tshift) - b./tau.*sin(tshift);
		else
			gammas.f.r =  ctf.*g./tau.^2;
			gammas.t.r = -g.*cft;
			%gammas.f.r =  g./tau.^2;
			%gammas.t.r = -g;
		end
	end
end

%% imaginary
if strcmp(only,'') || strcmp(only,'real')
	% recall that the imaginary part of Gamma affects the real part of K.
	if ~branch_flag
		[g,b,bc,tau,tshift,~,bsh] = branchparts_unpack(Sb,ids.gr > 0);
		if ids.Gami == 1
			% simplification means setting tau = 1
			tau = ones(size(tau));
		end
		if ids.gr > 0
			g = zeros(size(g));
		end
	else
		[g,b,bc,tau,tshift,~,bsh] = branchparts_unpack(Sb,ids.g);
		if ids.g > 0
			g = zeros(size(g));
		end
		if ids.simp == 1
			% simplification means setting tau = 1
			tau = ones(size(tau));
		end
	end
	b0 = bc/2;
	
	if ids.AB == 1
		gammas.f.m =  b./tau.*(T*ctmp).*cos(tshift) + g./tau.*(T*ctmp).*sin(tshift);
		gammas.t.m = -b./tau.*(F*ctmp).*cos(tshift) + g./tau.*(F*ctmp).*sin(tshift);
	else
		gammas.f.m =  b./tau.*(T*ctmp);
		gammas.t.m = -b.*(F*ctmp);
	end
end
