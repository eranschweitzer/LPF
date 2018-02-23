function S = zhigangflow(vars, F, T, E, Sb,varargin)

pq = varargin_parse(varargin, 'pq', '');

S = struct();
[g,b,bc,tau,tshift,~,~] = branchparts_unpack(Sb,0);

x = E*vars.u - log(tau) + 1i*(E*vars.theta - tshift);

ystar = (g - 1i*b);

S.f =  ystar.*conj(x) + 0.5*ystar.*x.*conj(x);
S.t = -ystar.*conj(x) + 0.5*ystar.*x.*conj(x); 

if strcmp(pq,'real')
	S.f = real(S.f);
	S.t = real(S.t);
elseif strcmp(pq, 'imag')
	S.f = imag(S.f) ;%- (bc/2).*(1 + F*vars.u).^2;
	S.t = imag(S.t) ;%- (bc/2).*(1 + T*vars.u).^2;
end
