function [matI,matU] = pv_ref_mats(bidx,N,u0,ids)

matI = struct(); matU = struct();
matI.pv  = sparse(1:N.pv,find(bidx.pv),1, N.pv, N.t);
matI.ref  = sparse(1:N.ref,find(bidx.ref),1, N.ref, N.t);

if nargin == 3
    ids = struct();
end
if (nargin==3) || (length(fieldnames(ids)) == 13)
    matU.pv  = sparse(find(bidx.pv),1:N.pv, exp(-u0(bidx.pv)), N.t, N.pv);
    matU.ref.q  = sparse(find(bidx.ref),1:N.ref, exp(-u0(bidx.ref)), N.t, N.ref);
    matU.ref.p  = sparse(find(bidx.ref),1:N.ref, -exp(-u0(bidx.ref)), N.t, N.ref);
else
%     if ids.Q == 2
%         aq = 2;
%     else
%         aq = 1;
%     end
% 		ap = 1;
%     if ids.c == 2
%         %c=2 neglects all initial voltage information therefore the uhat term is needed
% 				%note that since u0 is NOT 0 here for PV busses necesseraily, therefore it plays the role of u/uhat
%         auq = aq; 
%         aq   = 0;
%         ap   = 0;
%         matU.bu.pv = u0(bidx.pv);
%         matU.bu.ref= u0(bidx.ref);
%     else
%         auq = 0;
%     end
%     if (ids.c ~= 2) || (ids.gr == 2)
%         aup = 0;
%     else
%         aup = 1;
%     end
%     matU.pv  = sparse(find(bidx.pv),1:N.pv,   (1 - auq*u0(bidx.pv)).*exp(-aq*u0(bidx.pv)), N.t, N.pv);
%     matU.ref = sparse(find(bidx.ref),1:N.ref, (1- aup*u0(bidx.ref)).*exp(-ap*u0(bidx.ref)), N.t, N.ref);
    
    
    matU.pv     = sparse(find(bidx.pv), 1:N.pv,    exp(-ids.Q.*u0(bidx.pv)), N.t, N.pv);
    matU.ref.q  = sparse(find(bidx.ref),1:N.ref,   exp(-ids.Q.*u0(bidx.ref)), N.t, N.ref);
    matU.ref.p  = sparse(find(bidx.ref),1:N.ref,  -exp(-u0(bidx.ref)), N.t, N.ref);
    matU.bu.pv = u0(bidx.pv);
    matU.bu.ref= u0(bidx.ref);
end
