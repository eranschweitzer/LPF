function out = nearestu0(F,T,Sb,N,bidx)

% Adj = (F + T)'*(F+T); % adjacency matrix
% pqnodes = find(bidx.pq);
% pvnodes = find(bidx.pv);
% refnodes = find(bidx.ref);
% candidates = cell(N.pq,1);
% for ii = 1:N.pq
%     x = sparse(pqnodes(ii),1,1,N.t,1);
%     while ~any(x & (bidx.pv | bidx.ref))
%         x = x | Adj'*x; %BFS
%     end
%     tmp = find(x & (bidx.pv | bidx.ref));
%     
%     candidates{ii} = struct('pv',find(ismember(pvnodes,tmp)), 'ref',find(ismember(refnodes,tmp)) );
% end

Ytr = myMakeYbus(F,T,Sb,'tronly',true);
%order as pq buses then pv buses
Y   = [Ytr(bidx.pq, bidx.pq), Ytr(bidx.pq, bidx.pv);
       Ytr(bidx.pv, bidx.pq), Ytr(bidx.pv, bidx.pv)];
if N.t < 1e4
    Z   = Y\speye(N.pv + N.pq);
    zdiag = diag(Z);
    Z = abs(ones(N.pv+N.pq,1)*zdiag.' + zdiag*ones(1,N.pv+N.pq) - Z - Z.');
%     out = zeros(N.pq,1);
%     for ii = 1:N.pq
%         if ii == 1222
%             pqnodes(ii)
%         end
%         if ~isempty(candidates{ii}.pv)
%             %pick nearest electrical pv bus from nearest topological pv
%             %buses
%             [zpv,idxpv] = min(Z(ii,N.pq+candidates{ii}.pv));
%         else
%             zpv = inf;
%         end
%         if ~isempty(candidates{ii}.ref) && (abs(zdiag(ii)) < zpv)
%             out(ii) = 0;
%         else
%             out(ii) = candidates{ii}.pv(idxpv);
%         end
%     end
    [zpv,out] = min(Z(1:N.pq,N.pq+1:end),[],2);
    out(abs(zdiag(1:N.pq)) < zpv) = 0;
else
    [l,u] = lu(Y);

    out = zeros(N.pq,1);
    for k = 1:N.pq
        rhs  = sparse( [k*ones(N.pv,1); (N.pq+1:N.pq+N.pv)'], [1:N.pv,1:N.pv]', ...
            [ones(N.pv,1);-ones(N.pv,1)], N.pv+N.pq,N.pv);
        Ztmp = u\(l\rhs);
        Zth  = Ztmp(k,:) - diag(Ztmp(N.pq+1:N.pq+N.pv,1:N.pv)).';
        [zpv,pvidx] = min(abs(Zth));
        Zref = u\(l\sparse(k,1,1,N.pv+N.pq,1));
        if abs(Zref(k)) > zpv
            %%% closest is pv bus number pvidx
            %%% if Zref <= zpv then the closest bus is the reference which will
            %%% be indicated by a value of 0 in the output.
            out(k) = pvidx;
        end
    end
end