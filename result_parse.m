function S = result_parse(x,u0,N,varargin)

uvar = varargin_parse(varargin, 'uvar', 'uhat');

S = struct();
%parse result
S.theta = x(1:N.t);
S.(uvar)= x(N.t+1:2*N.t);
% S.u     = vlimit(S.uhat + u0);
% S.v     = exp(S.uhat + u0);
S.Qpv   = x(2*N.t+1:2*N.t+N.pv);
S.Pref  = x(2*N.t+N.pv+1:2*N.t+N.pv+N.ref);
S.Qref  = x(end-N.ref+1:end);
S.u0    = u0;
% if (nargin == 4) && isfield(ids,'c')
% 	if ids.c == 2
% 		%%%% if initial values were neglected in the formulation
% 		%%%% Note: in this situation uhat was forced to be u0 for PV and ref buses
% 		%%% thererefore for u = u0 + uhat to be correct u0 = 0
% 		S.u0 = zeros(size(u0));
% 	else
% 		S.u0 = u0;
% 	end
% else
% 	S.u0    = u0;
% end
