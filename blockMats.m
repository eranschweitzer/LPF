function [A,b,Ap] = blockMats(opmats,Sb,Sp,u0,tref,ids)
v2struct(opmats); %move fields to variables
ids = struct('Q',2, ...
    'grpq', 0, 'gipq', 0, 'grpv', 0, 'gipv', 0,...
    'brpq', 0, 'bipq', 0, 'brpv', 0, 'bipv', 0,...
    'trpq', 0, 'tipq', 0, 'trpv', 0, 'tipv', 0,...
    'urpq', 0, 'uipq', 0, 'urpv', 0, 'uipv', 0);
b = cell(4,1);
A = cell(4,4);
Ap= cell(4,1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% b vector components %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% real b PQ
yw = reweightY(opmats,Sb,Sp,u0,'pqopt',1, 'gs', ids.grpq > 0, 'rsimp', ids.brpq==1);
% b{1} = real(conj(yw.S(bidx.pq)) - yw.taupq - (yw.ysh(bidx.pq) + yw.yff(bidx.pq) + yw.ytt(bidx.pq) + ...
%             cnx.F.pq'*(yw.yft + 1i*tref*(yw.yft.*lidx.T.ref)) + ...
%             cnx.T.pq'*(yw.ytf.*(1 - lidx.F.pq.*yw.lntau) + 1i*tref*(yw.ytf.*lidx.F.ref))));
b{1} = real(conj(yw.S(bidx.pq)) + yw.taupq - (yw.ysh(bidx.pq) + yw.yff(bidx.pq) + yw.ytt(bidx.pq) + ...
            cnx.F.pq'*((1 + 1i*yw.delta + 1i*tref*lidx.T.ref).*yw.yft) + ...
            cnx.T.pq'*((1 - (yw.lntau + 1i*yw.delta) + 1i*tref*lidx.F.ref).*yw.ytf)...
            ));
        
%% imag b PQ
yw = reweightY(opmats,Sb,Sp,u0,'pqopt',ids.Q, 'gs', ids.gipq > 0, 'isimp', ids.bipq==1);
% b{2} = imag(conj(yw.S(bidx.pq)) - (yw.ysh(bidx.pq) + yw.yff(bidx.pq) + yw.ytt(bidx.pq) + ...
%             cnx.F.pq'*(yw.yft.*(1 + yw.lntau) + 1i*tref*(yw.yft.*lidx.T.ref)) + ...
%             cnx.T.pq'*(yw.ytf.*(1 - yw.lntau) + 1i*tref*(yw.ytf.*lidx.F.ref))));
b{2} = imag(conj(yw.S(bidx.pq)) - (yw.ysh(bidx.pq) + yw.yff(bidx.pq) + yw.ytt(bidx.pq) + ...
            cnx.F.pq'*((1 + (yw.lntau + 1i*yw.delta) + 1i*tref*lidx.T.ref).*yw.yft) + ...
            cnx.T.pq'*((1 - (yw.lntau + 1i*yw.delta) + 1i*tref*lidx.F.ref).*yw.ytf)...
            ));
%% real b PV
yw = reweightY(opmats,Sb,Sp,u0,'pqopt', 1, 'gs', ids.grpv > 0, 'rsimp', ids.brpv==1);
% b{3} = real(conj(yw.S(bidx.pv)) - (yw.ysh(bidx.pv) + yw.yff(bidx.pv) + yw.ytt(bidx.pv) + ...
%             cnx.F.pv'*(yw.yft + 1i*tref*(yw.yft.*lidx.T.ref)) + ...
%             cnx.T.pv'*(yw.ytf.*(1 - lidx.F.pq.*yw.lntau) + 1i*tref*(yw.ytf.*lidx.F.ref))));
% b{3} = real(conj(yw.S(bidx.pv)) - (yw.ysh(bidx.pv) + yw.yff(bidx.pv) + yw.ytt(bidx.pv) + ...
%             cnx.F.pv'*((1 + (yw.lntau + 1i*yw.delta) + 1i*tref*lidx.T.ref).*yw.yft) + ...
%             cnx.T.pv'*((1 - (yw.lntau + 1i*yw.delta) + 1i*tref*lidx.F.ref).*yw.ytf)...
%             ));
b{3} = real(conj(yw.S(bidx.pv)) + yw.taupv - (yw.ysh(bidx.pv) + yw.yff(bidx.pv) + yw.ytt(bidx.pv) + ...
            cnx.F.pv'*((1 + 1i*yw.delta + 1i*tref*lidx.T.ref).*yw.yft) + ...
            cnx.T.pv'*((1 - (yw.lntau + 1i*yw.delta) + 1i*tref*lidx.F.ref).*yw.ytf)...
            ));       
%% imag b PV
yw = reweightY(opmats,Sb,Sp,u0, 'pqopt', 2, 'gs', ids.gipv > 0, 'isimp', ids.bipv==1);
% b{4} = imag(conj(yw.S(bidx.pv)) - (yw.ysh(bidx.pv) + yw.yff(bidx.pv) + yw.ytt(bidx.pv) + ...
%             cnx.F.pv'*(yw.yft + 1i*tref*(yw.yft.*lidx.T.ref)) + ...
%             cnx.T.pv'*(yw.ytf.*(1 - lidx.F.pq.*yw.lntau) + 1i*tref*(yw.ytf.*lidx.F.ref))));
b{4} = imag(conj(yw.S(bidx.pv)) - (yw.ysh(bidx.pv) + yw.yff(bidx.pv) + yw.ytt(bidx.pv) + ...
            cnx.F.pv'*((1 + (yw.lntau + 1i*yw.delta) + 1i*tref*lidx.T.ref).*yw.yft) + ...
            cnx.T.pv'*((1 - (yw.lntau + 1i*yw.delta) + 1i*tref*lidx.F.ref).*yw.ytf)...
            ));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Real A PQ components %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% real A PQ theta (PQ and PV)
yw = reweightY(opmats,Sb,Sp,u0,'pqopt',1, 'gs', ids.grpq > 0, 'rsimp', ids.trpq==1);
A{1,1} = real(1i*( (-cnx.F.pq'*diags(yw.yft) + cnx.T.pq'*diags(yw.ytf.*lidx.F.pq))*cnx.F.pq +...
                (+cnx.F.pq'*diags(yw.yft.*lidx.T.pq) - cnx.T.pq'*diags(yw.ytf))*cnx.T.pq));
A{1,3} = real(1i*(cnx.F.pq'*diags(yw.yft.*lidx.T.pv)*cnx.T.pv + cnx.T.pq'*diags(yw.ytf.*lidx.F.pv)*cnx.F.pv));

%% real A PQ uhat
yw = reweightY(opmats,Sb,Sp,u0,'pqopt',1, 'gs', ids.grpq > 0, 'rsimp', ids.urpq==1);
A{1,2} = real(diags(conj(yw.S(bidx.pq)) + yw.ysh(bidx.pq) + yw.yff(bidx.pq) + yw.ytt(bidx.pq)) + ...
                    cnx.F.pq'*diags(yw.yft.*lidx.T.pq)*cnx.T.pq + ...
                    cnx.T.pq'*diags(yw.ytf.*lidx.F.pq)*cnx.F.pq);

%% real A PQ Qg
A{1,4} = sparse(N.pq,N.pv);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Imaginary A PQ components %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% imag A PQ theta (PQ and PV)
yw = reweightY(opmats,Sb,Sp,u0,'pqopt',ids.Q, 'gs', ids.gipq > 0, 'isimp', ids.tipq==1);
A{2,1} = imag(1i*( (-cnx.F.pq'*diags(yw.yft) + cnx.T.pq'*diags(yw.ytf.*lidx.F.pq))*cnx.F.pq +...
                (+cnx.F.pq'*diags(yw.yft.*lidx.T.pq) - cnx.T.pq'*diags(yw.ytf))*cnx.T.pq));
A{2,3} = imag(1i*(cnx.F.pq'*diags(yw.yft.*lidx.T.pv)*cnx.T.pv + cnx.T.pq'*diags(yw.ytf.*lidx.F.pv)*cnx.F.pv));

%% imag A PQ uhat
yw = reweightY(opmats,Sb,Sp,u0,'pqopt',ids.Q, 'gs', ids.gipq > 0, 'isimp', ids.uipq==1);
switch yw.pqopt
    case 1
        A{2,2} = imag(diags(conj(yw.S(bidx.pq)) + yw.ysh(bidx.pq) + yw.yff(bidx.pq) + yw.ytt(bidx.pq)) + ...
                    cnx.F.pq'*diags(yw.yft.*lidx.T.pq)*cnx.T.pq + ...
                    cnx.T.pq'*diags(yw.ytf.*lidx.F.pq)*cnx.F.pq);
    case 2
        A{2,2} = imag((-cnx.F.pq'*diags(yw.yft) + cnx.T.pq'*diags(yw.ytf.*lidx.F.pq))*cnx.F.pq +...
                (+cnx.F.pq'*diags(yw.yft.*lidx.T.pq) - cnx.T.pq'*diags(yw.ytf))*cnx.T.pq);
    otherwise
        error('yw.pqopt must be either 1 or 2')
end
%% imag A PQ Qg
A{2,4} = sparse(N.pq,N.pv); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Real A PV components %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% real A PV theta (PQ and PV)
yw = reweightY(opmats,Sb,Sp,u0, 'pqopt', 1, 'gs', ids.grpv > 0, 'rsimp', ids.trpv==1);
A{3,1} = real(1i*(cnx.F.pv'*diags(yw.yft.*lidx.T.pq)*cnx.T.pq + cnx.T.pv'*diags(yw.ytf.*lidx.F.pq)*cnx.F.pq));
A{3,3} = real(1i*( (-cnx.F.pv'*diags(yw.yft) + cnx.T.pv'*diags(yw.ytf.*lidx.F.pv))*cnx.F.pv +...
                (+cnx.F.pv'*diags(yw.yft.*lidx.T.pv) - cnx.T.pv'*diags(yw.ytf))*cnx.T.pv));
           
%% real A PV uhat
yw = reweightY(opmats,Sb,Sp,u0, 'pqopt', 1,'gs', ids.grpv > 0, 'rsimp', ids.urpv==1);
A{3,2} = real(cnx.F.pv'*diags(yw.yft.*lidx.T.pq)*cnx.T.pq + cnx.T.pv'*diags(yw.ytf.*lidx.F.pq)*cnx.F.pq);

%% real A PV Qg
A{3,4} = sparse(N.pv,N.pv);
% A{3,4} = sparse(1:N.pv,1:N.pv, exp(-u0(bidx.pv)),N.pv,N.pv);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Imaginary A PV components %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% imag A PV theta (PQ and PV)
yw = reweightY(opmats,Sb,Sp,u0, 'pqopt', 2, 'gs', ids.gipv > 0, 'isimp', ids.tipv==1,'vtau',false);
A{4,1} = imag(1i*(cnx.F.pv'*diags(yw.yft.*lidx.T.pq)*cnx.T.pq + cnx.T.pv'*diags(yw.ytf.*lidx.F.pq)*cnx.F.pq));
A{4,3} = imag(1i*( (-cnx.F.pv'*diags(yw.yft) + cnx.T.pv'*diags(yw.ytf.*lidx.F.pv))*cnx.F.pv +...
                (+cnx.F.pv'*diags(yw.yft.*lidx.T.pv) - cnx.T.pv'*diags(yw.ytf))*cnx.T.pv));
           
%% imag A PV uhat
yw = reweightY(opmats,Sb,Sp,u0, 'pqopt', 2,'gs', ids.gipv > 0, 'rsimp', ids.uipv==1,'vtau',false);
A{4,2} = imag(cnx.F.pv'*diags(yw.yft.*lidx.T.pq)*cnx.T.pq + cnx.T.pv'*diags(yw.ytf.*lidx.F.pq)*cnx.F.pq);

%% imag A PV Qg
A{4,4} = speye(N.pv,N.pv);
% A{4,4} = sparse(1:N.pv,1:N.pv, exp(-2*u0(bidx.pv)),N.pv,N.pv);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% A Phi components %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% real
yw = reweightY(opmats,Sb,Sp,u0);
Ap{1,1} = real(-(cnx.F.pq'*diags(yw.yft) + cnx.T.pq'*diags(yw.ytf)));
Ap{3,1} = real(-(cnx.F.pv'*diags(yw.yft) + cnx.T.pv'*diags(yw.ytf)));

%% imag
yw = reweightY(opmats,Sb,Sp,u0,'pqopt', ids.Q);
Ap{2,1} = imag(-(cnx.F.pq'*diags(yw.yft) + cnx.T.pq'*diags(yw.ytf)));
Ap{4,1} = imag(-(cnx.F.pv'*diags(yw.yft) + cnx.T.pv'*diags(yw.ytf)));
